#!/usr/bin/env python3
"""Corrected partial-Spearman correlation of metabolite features with color
phenotype, jointly adjusting for Species, Library Plate, and C/SUP sample type.

This replaces scripts/phase2_correlation_analysis.py for confidence purposes.
Two bugs in that script are fixed here, not just the missing-Species branch
documented in docs/FEATURE_ANALYSIS.md's caveats section:

1. `spearman_partial_corr()` looped over covariates but reassigned (not
   accumulated) the residual on each iteration, so with >1 covariate only the
   LAST covariate in the list was actually regressed out. This never surfaced
   because that run only ever passed one covariate (Library Plate); it would
   have silently under-corrected the moment Species was added. Fixed here with
   a single joint OLS regression against all covariates at once (rank-space
   projection), which is also the statistically correct way to compute a
   partial correlation with multiple covariates.
2. Library Plate was `pd.factorize`d into an arbitrary integer and regressed
   as if continuous/ordinal (plate "2" is not "twice" plate "1"). Fixed here
   by one-hot encoding Plate (and Species, and sample_type) as categorical
   dummies.

Method: rank-transform each feature column and each phenotype column
(Spearman), jointly regress out the covariate design matrix (OLS projection,
vectorized across all features via one matrix multiply), then Pearson-
correlate the residuals. This is the standard "partial Spearman via rank
residualization" approach and is equivalent, for a single covariate, to what
scripts/phase2_correlation_analysis.py intended to do.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata
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
RNG_SEED = 20260811  # documented per statistical-conventions.md reproducibility rule

print("=" * 80)
print("02: CORRECTED PARTIAL CORRELATION (Species + Plate + SampleType jointly regressed out)")
print("=" * 80)

design = pd.read_csv(os.path.join(OUT_DIR, "sample_design.csv"))
features = pd.read_csv(os.path.join(OUT_DIR, "features_cleaned.csv.gz"))
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert len(design) == len(features), "design/features row mismatch"
n = len(design)
print(f"samples: {n}  features: {features.shape[1]}  phenotypes: {len(PHENOTYPE_COLS)}")

# --- build joint design matrix (one-hot, drop-first, plus intercept) ---
Z_parts = [np.ones((n, 1))]
covariate_names = ["intercept"]
for col, prefix in [("species", "species"), ("Library Plate", "plate"), ("sample_type", "stype")]:
    dummies = pd.get_dummies(design[col].astype(str), prefix=prefix, drop_first=True)
    Z_parts.append(dummies.values.astype(float))
    covariate_names.extend(dummies.columns.tolist())
Z = np.hstack(Z_parts)
rank_Z = np.linalg.matrix_rank(Z)
print(f"design matrix: {Z.shape}, rank {rank_Z} (covariates: species={design['species'].nunique()}, "
      f"plate={design['Library Plate'].nunique()}, sample_type={design['sample_type'].nunique()})")

# Projection (hat) matrix via pseudo-inverse -- handles any residual collinearity
# between species/plate/sample_type by minimum-norm solution rather than crashing.
Z_pinv = np.linalg.pinv(Z)


def residualize(M):
    """Project out Z from every column of M (n x k) at once."""
    beta = Z_pinv @ M
    return M - Z @ beta


feature_values = features.values.astype(float)
feature_ranks = np.apply_along_axis(rankdata, 0, feature_values)
# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
feature_resid = residualize(feature_ranks)

# drop constant features (rank-degenerate after residualization, e.g. std==0 upstream)
feature_std = feature_resid.std(axis=0)
valid_feature_cols = np.where(feature_std > 1e-10)[0]
n_dropped_constant = feature_values.shape[1] - len(valid_feature_cols)
print(f"dropped {n_dropped_constant} features with ~zero residual variance after correction")

print("\n[1/3] Computing corrected partial correlations...")
results = []
for phenotype in PHENOTYPE_COLS:
    y = design[phenotype].values.astype(float)
    y_rank = rankdata(y)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
    y_resid = residualize(y_rank.reshape(-1, 1)).ravel()
    y_std = y_resid.std()
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert y_std > 1e-10, f"{phenotype} residual has ~zero variance after correction"

    # vectorized Pearson r of residuals = feature_resid^T y_resid / (norms), then p-values via scipy per-feature
    fr = feature_resid[:, valid_feature_cols]
    fr_c = fr - fr.mean(axis=0)
    y_c = y_resid - y_resid.mean()
    num = fr_c.T @ y_c
    denom = np.sqrt((fr_c ** 2).sum(axis=0) * (y_c ** 2).sum())
    rho = num / denom

    # degrees of freedom used up by the design matrix reduce effective n for the t-test
    df = n - rank_Z - 1
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert df > 5, f"too few residual degrees of freedom: {df}"
    t_stat = rho * np.sqrt(df / np.clip(1 - rho ** 2, 1e-12, None))
    from scipy.stats import t as t_dist
    pval = 2 * t_dist.sf(np.abs(t_stat), df=df)

    for j, feat_idx in enumerate(valid_feature_cols):
        results.append({
            "phenotype": phenotype,
            "feature_index": int(feat_idx),
            "rho_corrected": float(rho[j]),
            "pval_corrected": float(pval[j]),
            "n_samples": int(n),
            "df": int(df),
        })
    print(f"  {phenotype}: done ({len(valid_feature_cols)} features)")

df_results = pd.DataFrame(results)
print(f"\n✓ total corrected correlations: {len(df_results):,}")

print("\n[2/3] Two-stage BH-FDR correction...")
stage1 = []
for phenotype in PHENOTYPE_COLS:
    sub = df_results[df_results["phenotype"] == phenotype].copy()
    reject, q, _, _ = multipletests(sub["pval_corrected"], method="fdr_bh", alpha=0.05)
    sub["reject_stage1"] = reject
    sub["q_value_stage1"] = q
    stage1.append(sub)
    print(f"  {phenotype}: {reject.sum()} at q<0.05 (of {len(sub)})")
df_results = pd.concat(stage1, ignore_index=True)
reject_g, q_g, _, _ = multipletests(df_results["q_value_stage1"], method="fdr_bh", alpha=0.05)
df_results["q_value_global"] = q_g
df_results["reject_global"] = reject_g


def assign_tier(row):
    a = abs(row["rho_corrected"])
    if a > 0.30 and row["q_value_stage1"] < 0.05:
        return "Tier1_High"
    if a > 0.25 and row["q_value_stage1"] < 0.05:
        return "Tier2_Medium"
    if a > 0.20 and row["q_value_stage1"] < 0.10:
        return "Tier3_Exploratory"
    return "Not_Significant"


df_results["tier"] = df_results.apply(assign_tier, axis=1)
print("\nTiered results (corrected):")
print(df_results["tier"].value_counts().to_string())

print("\n[3/3] Comparing against the original (uncorrected) Phase 2 tier1 hit list...")
orig_tier1 = pd.read_csv(os.path.join(ROOT, "analysis", "phase2_tier1_hits.csv.gz"))
orig_tier1_keys = set(zip(orig_tier1["phenotype"], orig_tier1["feature_index"]))
corrected_tier1 = df_results[df_results["tier"] == "Tier1_High"]
corrected_tier1_keys = set(zip(corrected_tier1["phenotype"], corrected_tier1["feature_index"]))
overlap = orig_tier1_keys & corrected_tier1_keys
print(f"  original Tier1 hits: {len(orig_tier1_keys):,}")
print(f"  corrected Tier1 hits: {len(corrected_tier1_keys):,}")
print(f"  overlap: {len(overlap):,} ({100 * len(overlap) / max(len(orig_tier1_keys), 1):.1f}% of original survive)")

df_results.to_csv(os.path.join(OUT_DIR, "corrected_all_correlations.csv.gz"), index=False, compression="gzip")
corrected_tier1.sort_values("q_value_stage1").to_csv(
    os.path.join(OUT_DIR, "corrected_tier1_hits.csv.gz"), index=False, compression="gzip"
)
print(f"\n✓ wrote corrected_all_correlations.csv.gz, corrected_tier1_hits.csv.gz")

summary = {
    "n_samples": int(n),
    "n_features_tested": int(len(valid_feature_cols)),
    "n_features_dropped_constant": int(n_dropped_constant),
    "covariates": covariate_names,
    "design_matrix_rank": int(rank_Z),
    "residual_df": int(n - rank_Z - 1),
    "tier1_count_corrected": int((df_results["tier"] == "Tier1_High").sum()),
    "tier2_count_corrected": int((df_results["tier"] == "Tier2_Medium").sum()),
    "tier3_count_corrected": int((df_results["tier"] == "Tier3_Exploratory").sum()),
    "tier1_count_original_uncorrected": int(len(orig_tier1_keys)),
    "tier1_overlap_original_vs_corrected": int(len(overlap)),
}
with open(os.path.join(OUT_DIR, "corrected_correlation_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))

register_value("tier1_count_corrected", summary["tier1_count_corrected"], provenance="outputs/corrected_correlation_summary.json")
register_value("tier1_count_original_uncorrected", summary["tier1_count_original_uncorrected"], provenance="analysis/phase2_tier1_hits.csv.gz")
register_value("tier1_overlap_original_vs_corrected", summary["tier1_overlap_original_vs_corrected"], provenance="outputs/corrected_correlation_summary.json")
