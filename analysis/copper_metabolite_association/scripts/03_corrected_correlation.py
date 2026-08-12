#!/usr/bin/env python3
"""Species(+Plate)-corrected partial-Spearman correlation of metabolite features
with copper AUC, one track at a time (cell, supernatant -- independent feature
sets, see 01_prepare_data.py). Reuses the joint rank-space OLS residualization
method validated in
analysis/phenotype_metabolite_association/scripts/02_corrected_correlation.py
(single joint regression against all covariates at once, one-hot dummies for
categoricals, not sequential/single-covariate residualization).

Per-track covariates (see 02_confound_check.py): species is a strong confound in
BOTH tracks (Kruskal-Wallis p=5.8e-10). Library Plate is also confounded in the
cell track (p=9.6e-5) but is CONSTANT (=1.0 for all 264 supernatant rows -- an
upstream data property, not a bug) in the supernatant track, so it is included as
a covariate for cell only.
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

TRACK_COVARIATES = {
    "cell": ["species", "Library Plate"],
    "supernatant": ["species"],  # Library Plate constant for this track, see docstring
}
# ANALYSIS_OK[threshold]: Tier thresholds reproduced from
# analysis/phenotype_metabolite_association/scripts/02_corrected_correlation.py
# (itself reproduced from the legacy phase2_correlation_analysis.py) for direct
# comparability across this project's analyses.
TIER1_RHO, TIER1_Q = 0.30, 0.05
TIER2_RHO, TIER2_Q = 0.25, 0.05
TIER3_RHO, TIER3_Q = 0.20, 0.10

print("=" * 80)
print("03: CORRECTED PARTIAL CORRELATION (species [+plate] regressed out, per track)")
print("=" * 80)


def build_design(design_df, covariate_cols):
    n = len(design_df)
    Z_parts = [np.ones((n, 1))]
    names = ["intercept"]
    for col in covariate_cols:
        dummies = pd.get_dummies(design_df[col].astype(str), drop_first=True)
        Z_parts.append(dummies.values.astype(float))
        names.extend([f"{col}={c}" for c in dummies.columns])
    Z = np.hstack(Z_parts)
    return Z, np.linalg.matrix_rank(Z), np.linalg.pinv(Z), names


def residualize(M, Z, Z_pinv):
    return M - Z @ (Z_pinv @ M)


def assign_tier(rho, q):
    a = abs(rho)
    if a > TIER1_RHO and q < TIER1_Q:
        return "Tier1_High"
    if a > TIER2_RHO and q < TIER2_Q:
        return "Tier2_Medium"
    if a > TIER3_RHO and q < TIER3_Q:
        return "Tier3_Exploratory"
    return "Not_Significant"


track_summaries = {}
for track, covariate_cols in TRACK_COVARIATES.items():
    print("\n" + "-" * 80)
    print(f"TRACK: {track}  (covariates: {covariate_cols})")
    print("-" * 80)
    out_dir = os.path.join(OUT_DIR, track)

    design = pd.read_csv(os.path.join(out_dir, "sample_design.csv"))
    features = pd.read_csv(os.path.join(out_dir, "features_cleaned.csv.gz"))
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert len(design) == len(features), f"[{track}] design/features row mismatch"
    n = len(design)
    print(f"samples: {n}  features: {features.shape[1]}")

    Z, rank_Z, Z_pinv, covariate_names = build_design(design, covariate_cols)
    df = n - rank_Z - 1
    print(f"design matrix: {Z.shape}, rank {rank_Z}, residual df {df}")
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert df > 5, f"[{track}] too few residual degrees of freedom: {df}"

    feature_values = features.values.astype(float)
    feature_ranks = np.apply_along_axis(rankdata, 0, feature_values)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against species
    # (+plate for the cell track); method reused from
    # phenotype_metabolite_association/scripts/02_corrected_correlation.py.
    feature_resid = residualize(feature_ranks, Z, Z_pinv)

    feature_std = feature_resid.std(axis=0)
    valid_cols = np.where(feature_std > 1e-10)[0]
    n_dropped_constant = feature_values.shape[1] - len(valid_cols)
    print(f"dropped {n_dropped_constant} features with ~zero residual variance after correction")

    y = design["mean_auc_rate"].values.astype(float)
    y_rank = rankdata(y)
    # ANALYSIS_OK[batch-correction]: same rank-space residualization as the feature matrix above, applied to the phenotype.
    y_resid = residualize(y_rank.reshape(-1, 1), Z, Z_pinv).ravel()
    y_std = y_resid.std()
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert y_std > 1e-10, f"[{track}] copper AUC residual has ~zero variance after correction"

    fr = feature_resid[:, valid_cols]
    fr_c = fr - fr.mean(axis=0)
    y_c = y_resid - y_resid.mean()
    num = fr_c.T @ y_c
    denom = np.sqrt((fr_c ** 2).sum(axis=0) * (y_c ** 2).sum())
    rho = num / denom
    t_stat = rho * np.sqrt(df / np.clip(1 - rho ** 2, 1e-12, None))
    pval = 2 * t_dist.sf(np.abs(t_stat), df=df)

    results = pd.DataFrame({
        "feature_index": valid_cols.astype(int),
        "rho_corrected": rho,
        "pval_corrected": pval,
        "n_samples": n,
        "df": df,
    })

    reject, q, _, _ = multipletests(results["pval_corrected"], method="fdr_bh", alpha=0.05)
    results["q_value"] = q
    results["reject_fdr05"] = reject
    results["tier"] = [assign_tier(r, qq) for r, qq in zip(results["rho_corrected"], results["q_value"])]

    print("\nTiered results:")
    print(results["tier"].value_counts().to_string())
    print(f"\nmax |rho_corrected|: {results['rho_corrected'].abs().max():.4f}")
    print("\ntop 10 by pval_corrected:")
    print(results.sort_values("pval_corrected").head(10).to_string(index=False))

    results = results.sort_values("pval_corrected").reset_index(drop=True)
    results.to_csv(os.path.join(out_dir, "corrected_all_correlations.csv.gz"), index=False, compression="gzip")
    tier1 = results[results["tier"] == "Tier1_High"]
    tier1.to_csv(os.path.join(out_dir, "corrected_tier1_hits.csv.gz"), index=False, compression="gzip")
    print(f"\n✓ wrote {out_dir}/corrected_all_correlations.csv.gz, corrected_tier1_hits.csv.gz")

    track_summary = {
        "n_samples": int(n),
        "n_features_tested": int(len(valid_cols)),
        "n_features_dropped_constant": int(n_dropped_constant),
        "covariates": covariate_names,
        "design_matrix_rank": int(rank_Z),
        "residual_df": int(df),
        "max_abs_rho_corrected": float(results["rho_corrected"].abs().max()),
        "tier1_count": int((results["tier"] == "Tier1_High").sum()),
        "tier2_count": int((results["tier"] == "Tier2_Medium").sum()),
        "tier3_count": int((results["tier"] == "Tier3_Exploratory").sum()),
    }
    track_summaries[track] = track_summary
    for k in ("tier1_count", "tier2_count", "tier3_count", "max_abs_rho_corrected"):
        register_value(f"{track}_{k}", track_summary[k], provenance=f"outputs/{track}/corrected_all_correlations.csv.gz")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
with open(os.path.join(OUT_DIR, "corrected_correlation_summary.json"), "w") as f:
    json.dump(track_summaries, f, indent=2)
print(json.dumps(track_summaries, indent=2))
