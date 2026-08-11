#!/usr/bin/env python3
"""Permutation null for the corrected partial correlations (02_corrected_correlation.py).

Two reasons this is needed even though nothing cleared FDR<0.05 in step 2:
1. The parametric t-based p-value from Pearson-on-residuals assumes each row is
   an independent observation. It is not: each strain contributes up to two rows
   (C_*, SUP_*) with correlated phenotype residuals. This inflates the effective
   n and anti-conservatively shrinks p-values. A permutation null that respects
   the strain-level pairing gives a calibrated p-value instead.
2. It directly answers "is there ANY real signal, however small" for (a) the six
   headline features named in docs/FEATURE_ANALYSIS.md, and (b) the strongest
   nominal candidates from the corrected scan, rather than relying on FDR alone.

Permutation design (strain-block, respects the C/SUP pairing):
  - Strains with both a C_* and SUP_* row: their (y_resid_C, y_resid_SUP) pair is
    permuted as a unit across paired strains (shuffling which strain's covariate
    context -- species/plate/sample_type -- gets assigned which pair of outcome
    residuals).
  - Strains with only a C_* or only a SUP_* row: permuted separately within their
    own singleton stratum (same logic, block size 1).
This preserves the within-strain correlation of the two phenotype residuals and
the marginal distribution of y_resid, while breaking any link to the metabolite
feature values -- the correct null for "does this feature predict phenotype
beyond species/plate/sample_type, given the repeated-measures structure."
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

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
N_PERM = 5000
RNG_SEED = 20260811
TOP_N_PER_PHENOTYPE = 25
# The six headline features named in docs/FEATURE_ANALYSIS.md
HEADLINE_FEATURES = [2755, 6926, 6188, 1560, 5740, 2308]

print("=" * 80)
print("03: PERMUTATION NULL (strain-block, respects C/SUP pairing)")
print("=" * 80)

design = pd.read_csv(os.path.join(OUT_DIR, "sample_design.csv"))
features = pd.read_csv(os.path.join(OUT_DIR, "features_cleaned.csv.gz"))
corrected = pd.read_csv(os.path.join(OUT_DIR, "corrected_all_correlations.csv.gz"))
n = len(design)
rng = np.random.default_rng(RNG_SEED)

# --- rebuild the same design matrix / residualization as step 2 ---
Z_parts = [np.ones((n, 1))]
for col, prefix in [("species", "species"), ("Library Plate", "plate"), ("sample_type", "stype")]:
    dummies = pd.get_dummies(design[col].astype(str), prefix=prefix, drop_first=True)
    Z_parts.append(dummies.values.astype(float))
Z = np.hstack(Z_parts)
Z_pinv = np.linalg.pinv(Z)


def residualize(M):
    return M - Z @ (Z_pinv @ M)


feature_values = features.values.astype(float)
feature_ranks = np.apply_along_axis(rankdata, 0, feature_values)
# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
feature_resid_all = residualize(feature_ranks)

# --- candidate set: headline features + top-N nominal hits per phenotype ---
candidates = []
for phenotype in PHENOTYPE_COLS:
    sub = corrected[corrected["phenotype"] == phenotype].sort_values("pval_corrected")
    for feat_idx in sub.head(TOP_N_PER_PHENOTYPE)["feature_index"]:
        candidates.append((phenotype, int(feat_idx), "top_nominal"))
    for feat_idx in HEADLINE_FEATURES:
        candidates.append((phenotype, int(feat_idx), "headline_from_FEATURE_ANALYSIS.md"))
candidates_df = pd.DataFrame(candidates, columns=["phenotype", "feature_index", "source"]).drop_duplicates(
    subset=["phenotype", "feature_index"]
)
print(f"candidate feature-phenotype pairs to permutation-test: {len(candidates_df)}")

# --- strain-block permutation index builder ---
strain_of = design["strain_id"].values
type_of = design["sample_type"].values
row_idx = np.arange(n)

# map strain -> {'C': row or None, 'SUP': row or None}
strain_rows = {}
for i, s in enumerate(strain_of):
    strain_rows.setdefault(s, {})[type_of[i]] = i

paired_strains = [s for s, d in strain_rows.items() if set(d) == {"C", "SUP"}]
c_only_strains = [s for s, d in strain_rows.items() if set(d) == {"C"}]
sup_only_strains = [s for s, d in strain_rows.items() if set(d) == {"SUP"}]
print(f"strain-block strata: paired={len(paired_strains)}, C-only={len(c_only_strains)}, SUP-only={len(sup_only_strains)}")
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert len(paired_strains) + len(c_only_strains) + len(sup_only_strains) == len(strain_rows)
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert 2 * len(paired_strains) + len(c_only_strains) + len(sup_only_strains) == n


def permuted_row_order(rng):
    """Return an array `perm` (length n) such that y_resid[perm] is one valid
    strain-block permutation: row i's phenotype residual is replaced by the
    residual from the strain assigned to it within its own pairing stratum."""
    perm = np.empty(n, dtype=int)
    order = rng.permutation(len(paired_strains))
    for target_s, source_i in zip(paired_strains, order):
        source_s = paired_strains[source_i]
        for t in ("C", "SUP"):
            perm[strain_rows[target_s][t]] = strain_rows[source_s][t]
    for group in (c_only_strains, sup_only_strains):
        if not group:
            continue
        order = rng.permutation(len(group))
        for target_s, source_i in zip(group, order):
            source_s = group[source_i]
            t = next(iter(strain_rows[target_s]))
            perm[strain_rows[target_s][t]] = strain_rows[source_s][t]
    return perm


perm_results = []
for phenotype in PHENOTYPE_COLS:
    y = design[phenotype].values.astype(float)
    y_rank = rankdata(y)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
    y_resid = residualize(y_rank.reshape(-1, 1)).ravel()

    feat_idxs = candidates_df.loc[candidates_df["phenotype"] == phenotype, "feature_index"].unique()
    fr = feature_resid_all[:, feat_idxs]
    fr_c = fr - fr.mean(axis=0)
    fr_norm = np.sqrt((fr_c ** 2).sum(axis=0))

    y_c = y_resid - y_resid.mean()
    obs_rho = (fr_c.T @ y_c) / (fr_norm * np.sqrt((y_c ** 2).sum()))

    null_rho = np.empty((N_PERM, len(feat_idxs)))
    for p in range(N_PERM):
        perm = permuted_row_order(rng)
        yp = y_resid[perm]
        yp_c = yp - yp.mean()
        null_rho[p, :] = (fr_c.T @ yp_c) / (fr_norm * np.sqrt((yp_c ** 2).sum()))

    perm_p = (np.sum(np.abs(null_rho) >= np.abs(obs_rho), axis=0) + 1) / (N_PERM + 1)

    for j, feat_idx in enumerate(feat_idxs):
        perm_results.append({
            "phenotype": phenotype,
            "feature_index": int(feat_idx),
            "rho_corrected": float(obs_rho[j]),
            "perm_pval": float(perm_p[j]),
            "null_rho_2p5": float(np.percentile(null_rho[:, j], 2.5)),
            "null_rho_97p5": float(np.percentile(null_rho[:, j], 97.5)),
        })
    print(f"  {phenotype}: permutation-tested {len(feat_idxs)} features x {N_PERM} perms")

# ANALYSIS_OK[join]: perm_results has exactly one row per (phenotype, feature_index)
# already deduplicated in candidates_df; validate enforces that.
perm_df = pd.DataFrame(perm_results).merge(
    candidates_df, on=["phenotype", "feature_index"], how="left", validate="one_to_one"
)
perm_df = perm_df.sort_values("perm_pval")
perm_df.to_csv(os.path.join(OUT_DIR, "permutation_null_results.csv"), index=False)

n_headline_significant = int(
    ((perm_df["source"] == "headline_from_FEATURE_ANALYSIS.md") & (perm_df["perm_pval"] < 0.05)).sum()
) if len(perm_df) else 0
n_candidates_significant = int((perm_df["perm_pval"] < 0.05).sum())

print(f"\n✓ wrote permutation_null_results.csv")
print(f"headline features (FEATURE_ANALYSIS.md) with perm_pval<0.05: {n_headline_significant}/{len(HEADLINE_FEATURES) * len(PHENOTYPE_COLS)} pairs tested")
print(f"any candidate (headline + top-nominal) with perm_pval<0.05: {n_candidates_significant}/{len(perm_df)}")
print("\nTop 15 by permutation p-value:")
print(perm_df.head(15).to_string(index=False))
print("\nHeadline features specifically:")
print(perm_df[perm_df["source"] == "headline_from_FEATURE_ANALYSIS.md"].to_string(index=False))

summary = {
    "n_permutations": N_PERM,
    "n_candidates_tested": int(len(perm_df)),
    "n_candidates_perm_significant_p05": n_candidates_significant,
    "n_headline_features_tested": int((perm_df["source"] == "headline_from_FEATURE_ANALYSIS.md").sum()),
    "n_headline_features_perm_significant_p05": n_headline_significant,
    "rng_seed": RNG_SEED,
}
with open(os.path.join(OUT_DIR, "permutation_null_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))

register_value("n_permutations", N_PERM, provenance="scripts/03_permutation_null.py")
register_value("n_candidates_perm_significant_p05", n_candidates_significant, provenance="outputs/permutation_null_results.csv")
register_value("n_headline_features_perm_significant_p05", n_headline_significant, provenance="outputs/permutation_null_results.csv")
