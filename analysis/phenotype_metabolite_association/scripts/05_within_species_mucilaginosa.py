#!/usr/bin/env python3
"""Within-species re-analysis restricted to R. mucilaginosa, the only species in
this dataset with enough strains (n=210) to look for a phenotype-metabolite
association that isn't just "these are different species." Species is dropped
from the covariate set (constant within this subset); Library Plate and C/SUP
sample_type are still regressed out. Uses the same rank-residualization +
strain-block permutation approach as 02/03, restricted to this one species.

Rationale: a real, single-locus/pathway-level genetic effect on pigmentation
should be detectable *within* a species with a genetically diverse strain
collection, without needing between-species differences to generate the
correlation. docs/PHASE3_STRATIFIED_ANALYSIS_SUMMARY.md already found 0%
significant within-species hits for R. mucilaginosa using the OLD (buggy,
single-covariate) method on a smaller top-200-feature discriminant subset; this
re-checks with the corrected method across the full 7,341-feature set, plus a
permutation null.
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
SPECIES = "Rhodotorula mucilaginosa"
N_PERM = 5000
RNG_SEED = 20260811
TOP_N_PER_PHENOTYPE = 25
HEADLINE_FEATURES = [2755, 6926, 6188, 1560, 5740, 2308]

print("=" * 80)
print(f"05: WITHIN-SPECIES RE-ANALYSIS ({SPECIES})")
print("=" * 80)

design = pd.read_csv(os.path.join(OUT_DIR, "sample_design.csv"))
features = pd.read_csv(os.path.join(OUT_DIR, "features_cleaned.csv.gz"))
mask = (design["species"] == SPECIES).values
design_sp = design.loc[mask].reset_index(drop=True)
features_sp = features.loc[mask].reset_index(drop=True)
n = len(design_sp)
n_strains = design_sp["strain_id"].nunique()
print(f"samples: {n}  strains: {n_strains}  features: {features_sp.shape[1]}")

print("\nPhenotype spread within this species (is there enough variance to find anything?):")
spread_rows = []
for col in PHENOTYPE_COLS:
    v = design_sp[col]
    row = {"phenotype": col, "mean": v.mean(), "sd": v.std(), "min": v.min(), "max": v.max(), "cv_abs": v.std() / abs(v.mean())}
    spread_rows.append(row)
    print(f"  {col:28} mean={row['mean']:7.3f}  sd={row['sd']:6.3f}  range=[{row['min']:7.3f}, {row['max']:7.3f}]  |CV|={row['cv_abs']:.3f}")
spread_df = pd.DataFrame(spread_rows)
spread_df.to_csv(os.path.join(OUT_DIR, "mucilaginosa_phenotype_spread.csv"), index=False)

# --- design matrix: Library Plate + sample_type only (species constant here) ---
Z_parts = [np.ones((n, 1))]
covariate_names = ["intercept"]
for col, prefix in [("Library Plate", "plate"), ("sample_type", "stype")]:
    dummies = pd.get_dummies(design_sp[col].astype(str), prefix=prefix, drop_first=True)
    Z_parts.append(dummies.values.astype(float))
    covariate_names.extend(dummies.columns.tolist())
Z = np.hstack(Z_parts)
rank_Z = np.linalg.matrix_rank(Z)
Z_pinv = np.linalg.pinv(Z)
print(f"\ndesign matrix: {Z.shape}, rank {rank_Z}, covariates: {covariate_names}")


def residualize(M):
    return M - Z @ (Z_pinv @ M)


fvals = features_sp.values.astype(float)
franks = np.apply_along_axis(rankdata, 0, fvals)
# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
fresid = residualize(franks)
fresid_c = fresid - fresid.mean(axis=0)
fnorm = np.sqrt((fresid_c ** 2).sum(axis=0))
valid_cols = np.where(fnorm > 1e-10)[0]
print(f"dropped {fvals.shape[1] - len(valid_cols)} near-constant features within this species subset")

df_resid = n - rank_Z - 1
print(f"\n[1/2] Computing within-species corrected partial correlations (df={df_resid})...")
results = []
for phenotype in PHENOTYPE_COLS:
    y = design_sp[phenotype].values.astype(float)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
    y_resid = residualize(rankdata(y).reshape(-1, 1)).ravel()
    y_c = y_resid - y_resid.mean()

    fr = fresid_c[:, valid_cols]
    fn = fnorm[valid_cols]
    rho = (fr.T @ y_c) / (fn * np.sqrt((y_c ** 2).sum()))
    t_stat = rho * np.sqrt(df_resid / np.clip(1 - rho ** 2, 1e-12, None))
    pval = 2 * t_dist.sf(np.abs(t_stat), df=df_resid)
    for j, feat_idx in enumerate(valid_cols):
        results.append({"phenotype": phenotype, "feature_index": int(feat_idx), "rho": float(rho[j]), "pval": float(pval[j])})
    print(f"  {phenotype}: done")

df_results = pd.DataFrame(results)
stage1 = []
for phenotype in PHENOTYPE_COLS:
    sub = df_results[df_results["phenotype"] == phenotype].copy()
    reject, q, _, _ = multipletests(sub["pval"], method="fdr_bh", alpha=0.05)
    sub["q_value_stage1"] = q
    sub["reject_stage1"] = reject
    stage1.append(sub)
    print(f"  {phenotype}: {reject.sum()} at q<0.05 (of {len(sub)})")
df_results = pd.concat(stage1, ignore_index=True)
reject_g, q_g, _, _ = multipletests(df_results["q_value_stage1"], method="fdr_bh", alpha=0.05)
df_results["q_value_global"] = q_g

df_results.to_csv(os.path.join(OUT_DIR, "mucilaginosa_all_correlations.csv.gz"), index=False, compression="gzip")
print(f"max |rho| overall: {df_results['rho'].abs().max():.3f}")
print("\nTop 10 by nominal p-value:")
print(df_results.sort_values("pval").head(10).to_string(index=False))

# --- strain-block permutation, same design as 03_permutation_null.py, single-species subset ---
print("\n[2/2] Strain-block permutation null (top-N per phenotype + headline features)...")
strain_of = design_sp["strain_id"].values
type_of = design_sp["sample_type"].values
strain_rows = {}
for i, s in enumerate(strain_of):
    strain_rows.setdefault(s, {})[type_of[i]] = i
paired = [s for s, d in strain_rows.items() if set(d) == {"C", "SUP"}]
c_only = [s for s, d in strain_rows.items() if set(d) == {"C"}]
sup_only = [s for s, d in strain_rows.items() if set(d) == {"SUP"}]
print(f"strain-block strata: paired={len(paired)}, C-only={len(c_only)}, SUP-only={len(sup_only)}")
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


candidates = []
for phenotype in PHENOTYPE_COLS:
    sub = df_results[df_results["phenotype"] == phenotype].sort_values("pval")
    for feat_idx in sub.head(TOP_N_PER_PHENOTYPE)["feature_index"]:
        candidates.append((phenotype, int(feat_idx), "top_nominal"))
    for feat_idx in HEADLINE_FEATURES:
        candidates.append((phenotype, int(feat_idx), "headline_from_FEATURE_ANALYSIS.md"))
cand_df = pd.DataFrame(candidates, columns=["phenotype", "feature_index", "source"]).drop_duplicates(
    subset=["phenotype", "feature_index"]
)

perm_rows = []
for phenotype in PHENOTYPE_COLS:
    y = design_sp[phenotype].values.astype(float)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
    y_resid = residualize(rankdata(y).reshape(-1, 1)).ravel()
    feat_idxs = cand_df.loc[cand_df["phenotype"] == phenotype, "feature_index"].unique()
    fr = fresid[:, feat_idxs]
    fr_c = fr - fr.mean(axis=0)
    fr_norm = np.sqrt((fr_c ** 2).sum(axis=0))
    y_c = y_resid - y_resid.mean()
    obs_rho = (fr_c.T @ y_c) / (fr_norm * np.sqrt((y_c ** 2).sum()))

    null_rho = np.empty((N_PERM, len(feat_idxs)))
    for p in range(N_PERM):
        perm = permuted_row_order()
        yp = y_resid[perm]
        yp_c = yp - yp.mean()
        null_rho[p, :] = (fr_c.T @ yp_c) / (fr_norm * np.sqrt((yp_c ** 2).sum()))
    perm_p = (np.sum(np.abs(null_rho) >= np.abs(obs_rho), axis=0) + 1) / (N_PERM + 1)
    for j, feat_idx in enumerate(feat_idxs):
        perm_rows.append({"phenotype": phenotype, "feature_index": int(feat_idx), "rho": float(obs_rho[j]), "perm_pval": float(perm_p[j])})
    print(f"  {phenotype}: permutation-tested {len(feat_idxs)} features x {N_PERM} perms")

perm_df = pd.DataFrame(perm_rows).merge(cand_df, on=["phenotype", "feature_index"], how="left", validate="one_to_one").sort_values("perm_pval")
perm_df.to_csv(os.path.join(OUT_DIR, "mucilaginosa_permutation_results.csv"), index=False)

n_sig = int((perm_df["perm_pval"] < 0.05).sum())
n_headline_sig = int(((perm_df["source"] == "headline_from_FEATURE_ANALYSIS.md") & (perm_df["perm_pval"] < 0.05)).sum())
print(f"\ncandidates with perm_pval<0.05: {n_sig}/{len(perm_df)}")
print(f"headline features with perm_pval<0.05: {n_headline_sig}")
print("\nTop 15 within-species candidates by permutation p:")
print(perm_df.head(15).to_string(index=False))

summary = {
    "species": SPECIES,
    "n_samples": int(n),
    "n_strains": int(n_strains),
    "phenotype_spread": spread_df.to_dict(orient="records"),
    "n_features_tested": int(len(valid_cols)),
    "tier1_fdr_hits": int((df_results["q_value_stage1"] < 0.05).sum()),
    "max_abs_rho": float(df_results["rho"].abs().max()),
    "n_permutations": N_PERM,
    "n_candidates_perm_significant_p05": n_sig,
    "n_headline_perm_significant_p05": n_headline_sig,
    "rng_seed": RNG_SEED,
}
with open(os.path.join(OUT_DIR, "mucilaginosa_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n" + json.dumps(summary, indent=2))

register_value("mucilaginosa_n_samples", int(n), provenance="outputs/mucilaginosa_summary.json")
register_value("mucilaginosa_n_strains", int(n_strains), provenance="outputs/mucilaginosa_summary.json")
register_value("mucilaginosa_tier1_fdr_hits", summary["tier1_fdr_hits"], provenance="outputs/mucilaginosa_summary.json")
register_value("mucilaginosa_n_candidates_perm_significant_p05", n_sig, provenance="outputs/mucilaginosa_summary.json")
