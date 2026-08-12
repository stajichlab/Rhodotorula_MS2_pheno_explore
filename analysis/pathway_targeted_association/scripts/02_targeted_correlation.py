#!/usr/bin/env python3
"""Corrected partial correlation + strain-block permutation test for the pigment-
pathway-matched features from 01_match_targets.py, using the exact same
confound-correction (Species + Library Plate + sample_type regressed out jointly) and
permutation design as analysis/phenotype_metabolite_association/scripts/02 and 03 --
but pulling raw peak areas directly from the aligned table by `raw_position`, since
2/3 matched features were excluded from features_cleaned.csv.gz by Phase 1's generic
QC filters (see 01_match_targets.py docstring).

Normalization: total-ion-signal per sample (sum of ALL 16,332 raw feature peak areas
for that sample), NOT the sum-of-the-7,341-QC-survivors normalization Phase 1 used --
that denominator is only defined for the pre-filtered subset and would be
inconsistent to reuse for features Phase 1 excluded. Using the full-table total is
the standard "total ion current" normalization and doesn't depend on which features
happened to pass an unrelated QC filter.

With only 3 candidate features x 3 phenotypes = 9 tests (vs. 22,023 in the untargeted
scan), the BH-FDR multiple-testing burden is negligible -- this is the entire point
of going targeted, and directly improves the detection floor found by
09_power_analysis.py.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as t_dist
from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "outputs")
PMA_OUT_DIR = os.path.join(ROOT, "analysis", "phenotype_metabolite_association", "outputs")

sys.path.insert(0, os.path.join(ROOT, "skills", "core", "scripts"))
try:
    from register_value import register_value
except ImportError:
    # ANALYSIS_OK[optional-dependency]: register_value is an internal mycelium
    # reporting helper, not a scientific dependency; numbers are still printed to
    # stdout/saved in outputs/*.json above, so a no-op fallback loses convenience
    # (report auto-fill) but not any result.
    def register_value(*args, **kwargs):
        pass

PHENOTYPE_COLS = ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
RAW_TABLE_PATH = os.path.join(ROOT, "input_data", "Rhodotorula_MS2_aligned_features_ms2.csv.gz")
N_PERM = 5000
RNG_SEED = 20260811

print("=" * 80)
print("02: TARGETED CORRECTED CORRELATION (pigment-pathway candidates only)")
print("=" * 80)

targets = pd.read_csv(os.path.join(OUT_DIR, "pathway_target_feature_list.csv"))
design = pd.read_csv(os.path.join(PMA_OUT_DIR, "sample_design.csv"))
n = len(design)
print(f"candidate features: {len(targets)}  samples: {n}")

# --- pull raw peak areas for the candidate features + compute per-sample total ion signal ---
raw = pd.read_csv(RAW_TABLE_PATH)
sample_cols = [c for c in raw.columns if "Peak area" in c]
sample_total = raw[sample_cols].sum(axis=0)  # total ion signal per sample, ALL 16,332 features

candidate_raw = raw.iloc[targets["raw_position"].values][sample_cols]
candidate_raw.index = targets["raw_position"].values

# map design's sample_id ("C_1") to raw table's column name ("C_1.mzML Peak area")
sample_col_of = {c.replace(".mzML Peak area", ""): c for c in sample_cols}
missing = set(design["sample_id"]) - set(sample_col_of)
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert not missing, f"{len(missing)} design sample_ids not found in raw table columns: {sorted(missing)[:5]}"
design_cols = [sample_col_of[s] for s in design["sample_id"]]

abundance = candidate_raw[design_cols].T  # rows=samples (design order), cols=raw_position
abundance.index = design["sample_id"].values
totals = sample_total[design_cols].values  # aligned to design order

normalized = abundance.div(totals, axis=0)
pseudocount = 0.1 * normalized[normalized > 0].median().median()
log_abundance = np.log2(normalized + pseudocount)
print(f"pseudocount: {pseudocount:.3e}")

detection_rate = (abundance > 0).mean(axis=0)
print("\nDetection rate (fraction of 550 samples with nonzero peak area):")
for raw_pos, rate in detection_rate.items():
    compound = targets.loc[targets["raw_position"] == raw_pos, "compound"].iloc[0]
    print(f"  raw_position={raw_pos} ({compound}): {rate:.1%}")

# --- same design matrix / residualization as 02_corrected_correlation.py ---
Z_parts = [np.ones((n, 1))]
covariate_names = ["intercept"]
for col in ("species", "Library Plate", "sample_type"):
    dummies = pd.get_dummies(design[col].astype(str), drop_first=True)
    Z_parts.append(dummies.values.astype(float))
    covariate_names.extend(dummies.columns.tolist())
Z = np.hstack(Z_parts)
rank_Z = np.linalg.matrix_rank(Z)
Z_pinv = np.linalg.pinv(Z)
df_resid = n - rank_Z - 1


def residualize(M):
    return M - Z @ (Z_pinv @ M)


feature_ranks = log_abundance.apply(rankdata, axis=0).values
# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type, identical method to phenotype_metabolite_association/scripts/02_corrected_correlation.py.
feature_resid = residualize(feature_ranks)

print(f"\n[1/2] Corrected partial correlation (df={df_resid})...")
results = []
for phenotype in PHENOTYPE_COLS:
    y = design[phenotype].values.astype(float)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type, identical method to phenotype_metabolite_association/scripts/02_corrected_correlation.py.
    y_resid = residualize(rankdata(y).reshape(-1, 1)).ravel()
    y_c = y_resid - y_resid.mean()

    fr = feature_resid - feature_resid.mean(axis=0)
    fn = np.sqrt((fr ** 2).sum(axis=0))
    rho = (fr.T @ y_c) / (fn * np.sqrt((y_c ** 2).sum()))
    t_stat = rho * np.sqrt(df_resid / np.clip(1 - rho ** 2, 1e-12, None))
    pval = 2 * t_dist.sf(np.abs(t_stat), df=df_resid)

    for j, raw_pos in enumerate(log_abundance.columns):
        compound = targets.loc[targets["raw_position"] == raw_pos, "compound"].iloc[0]
        results.append({
            "phenotype": phenotype, "raw_position": int(raw_pos), "compound": compound,
            "rho": float(rho[j]), "pval": float(pval[j]),
        })

df_results = pd.DataFrame(results)
reject, q, _, _ = multipletests(df_results["pval"], method="fdr_bh", alpha=0.05)
df_results["q_value"] = q
df_results["significant_q05"] = reject
print(df_results.sort_values("pval").to_string(index=False))

# --- strain-block permutation, identical design to 03_permutation_null.py ---
print(f"\n[2/2] Strain-block permutation null ({N_PERM} permutations per phenotype)...")
strain_of = design["strain_id"].values
type_of = design["sample_type"].values
strain_rows = {}
for i, s in enumerate(strain_of):
    strain_rows.setdefault(s, {})[type_of[i]] = i
paired = [s for s, d in strain_rows.items() if set(d) == {"C", "SUP"}]
c_only = [s for s, d in strain_rows.items() if set(d) == {"C"}]
sup_only = [s for s, d in strain_rows.items() if set(d) == {"SUP"}]
rng = np.random.default_rng(RNG_SEED)


def permuted_row_order():
    perm = np.empty(n, dtype=int)
    order = rng.permutation(len(paired))
    for target_s, src_i in zip(paired, order):
        src_s = paired[src_i]
        for t in ("C", "SUP"):
            perm[strain_rows[target_s][t]] = strain_rows[src_s][t]
    for group in (c_only, sup_only):
        if not group:
            continue
        order = rng.permutation(len(group))
        for target_s, src_i in zip(group, order):
            src_s = group[src_i]
            t = next(iter(strain_rows[target_s]))
            perm[strain_rows[target_s][t]] = strain_rows[src_s][t]
    return perm


perm_rows = []
fr_c = feature_resid - feature_resid.mean(axis=0)
fr_norm = np.sqrt((fr_c ** 2).sum(axis=0))
for phenotype in PHENOTYPE_COLS:
    y = design[phenotype].values.astype(float)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type, identical method to phenotype_metabolite_association/scripts/02_corrected_correlation.py.
    y_resid = residualize(rankdata(y).reshape(-1, 1)).ravel()
    obs_rho = df_results.loc[df_results["phenotype"] == phenotype, "rho"].values

    null_rho = np.empty((N_PERM, len(log_abundance.columns)))
    for p in range(N_PERM):
        perm = permuted_row_order()
        yp = y_resid[perm]
        yp_c = yp - yp.mean()
        null_rho[p, :] = (fr_c.T @ yp_c) / (fr_norm * np.sqrt((yp_c ** 2).sum()))

    perm_p = (np.sum(np.abs(null_rho) >= np.abs(obs_rho), axis=0) + 1) / (N_PERM + 1)
    for j, raw_pos in enumerate(log_abundance.columns):
        compound = targets.loc[targets["raw_position"] == raw_pos, "compound"].iloc[0]
        perm_rows.append({"phenotype": phenotype, "raw_position": int(raw_pos), "compound": compound,
                           "rho": float(obs_rho[j]), "perm_pval": float(perm_p[j])})

perm_df = pd.DataFrame(perm_rows).sort_values("perm_pval")
print(perm_df.to_string(index=False))

df_results.to_csv(os.path.join(OUT_DIR, "targeted_correlation_results.csv"), index=False)
perm_df.to_csv(os.path.join(OUT_DIR, "targeted_permutation_results.csv"), index=False)

n_sig_fdr = int(df_results["significant_q05"].sum())
n_sig_perm = int((perm_df["perm_pval"] < 0.05).sum())
summary = {
    "n_candidate_features": int(len(targets)),
    "n_tests": int(len(df_results)),
    "n_significant_q05": n_sig_fdr,
    "n_significant_permutation_p05": n_sig_perm,
    "n_permutations": N_PERM,
    "detection_rates": {str(k): float(v) for k, v in detection_rate.items()},
    "rng_seed": RNG_SEED,
}
with open(os.path.join(OUT_DIR, "targeted_correlation_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n" + json.dumps(summary, indent=2))

register_value("targeted_n_significant_q05", n_sig_fdr, provenance="outputs/targeted_correlation_summary.json")
register_value("targeted_n_significant_permutation_p05", n_sig_perm, provenance="outputs/targeted_correlation_summary.json")
