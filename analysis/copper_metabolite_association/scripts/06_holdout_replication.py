#!/usr/bin/env python3
"""Species-stratified, strain-level holdout replication check, per track. Simpler
than analysis/phenotype_metabolite_association/scripts/04_replication_holdout.py's
version: each row here is already one strain (no C/SUP pairing within a track), so
splitting by row = splitting by strain directly.

Species with too few strains to split reliably (<5) are kept entirely in train and
are not part of the held-out test set (documented, not silently dropped).

Two checks, same as the color-phenotype analysis:
  1. Global calibration: Spearman rho of rho_train vs rho_test across ALL features.
  2. Hit replication: do the top-N nominal TRAIN hits replicate (same sign AND
     nominal p<0.05) in TEST?
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as t_dist

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
    "supernatant": ["species"],
}
RNG_SEED = 20260811
TEST_FRAC = 0.20
MIN_STRAINS_FOR_SPLIT = 5
TOP_N = 25

print("=" * 80)
print("06: SPECIES-STRATIFIED, STRAIN-LEVEL HOLDOUT REPLICATION (per track)")
print("=" * 80)


def fit_corrected(design, features, mask, covariate_cols):
    sub_design = design.loc[mask].reset_index(drop=True)
    sub_features = features.loc[mask].reset_index(drop=True)
    nn = len(sub_design)

    Z_parts = [np.ones((nn, 1))]
    for col in covariate_cols:
        dummies = pd.get_dummies(sub_design[col].astype(str), drop_first=True)
        Z_parts.append(dummies.values.astype(float))
    Z = np.hstack(Z_parts)
    rank_Z = np.linalg.matrix_rank(Z)
    Z_pinv = np.linalg.pinv(Z)

    def residualize(M):
        return M - Z @ (Z_pinv @ M)

    fvals = sub_features.values.astype(float)
    franks = np.apply_along_axis(rankdata, 0, fvals)
    # ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py, refit on this split only.
    fresid = residualize(franks)
    fresid_c = fresid - fresid.mean(axis=0)
    fnorm = np.sqrt((fresid_c ** 2).sum(axis=0))
    fnorm[fnorm < 1e-10] = np.nan  # constant features -> NaN corr, filtered later

    df = nn - rank_Z - 1
    y = sub_design["mean_auc_rate"].values.astype(float)
    # ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py, refit on this split only.
    y_resid = residualize(rankdata(y).reshape(-1, 1)).ravel()
    y_c = y_resid - y_resid.mean()
    rho = (fresid_c.T @ y_c) / (fnorm * np.sqrt((y_c ** 2).sum()))
    t_stat = rho * np.sqrt(df / np.clip(1 - rho ** 2, 1e-12, None))
    pval = 2 * t_dist.sf(np.abs(t_stat), df=df)
    return pd.DataFrame({"feature_index": np.arange(fvals.shape[1]), "rho": rho, "pval": pval}), df


overall_summary = {}
for track, covariate_cols in TRACK_COVARIATES.items():
    print("\n" + "-" * 80)
    print(f"TRACK: {track}")
    print("-" * 80)
    out_dir = os.path.join(OUT_DIR, track)
    rng = np.random.default_rng(RNG_SEED)

    design = pd.read_csv(os.path.join(out_dir, "sample_design.csv"))
    features = pd.read_csv(os.path.join(out_dir, "features_cleaned.csv.gz"))

    species_counts = design["species"].value_counts()
    splittable_species = species_counts[species_counts >= MIN_STRAINS_FOR_SPLIT].index
    print(f"species with >={MIN_STRAINS_FOR_SPLIT} strains (eligible for holdout split): "
          f"{len(splittable_species)}/{len(species_counts)}")

    test_idx = set()
    for sp in splittable_species:
        rows_sp = design.index[design["species"] == sp].tolist()
        n_test = max(1, round(len(rows_sp) * TEST_FRAC))
        chosen = rng.choice(rows_sp, size=n_test, replace=False)
        test_idx.update(chosen.tolist())

    is_test = design.index.isin(test_idx)
    n_test_rows, n_train_rows = int(is_test.sum()), int((~is_test).sum())
    print(f"train rows: {n_train_rows}  test rows: {n_test_rows}")
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert n_test_rows > 30, f"[{track}] held-out set too small to be informative"

    train_fit, train_df = fit_corrected(design, features, ~is_test, covariate_cols)
    test_fit, test_df = fit_corrected(design, features, is_test, covariate_cols)
    print(f"train residual df={train_df}, test residual df={test_df}")

    tr = train_fit.rename(columns={"rho": "rho_train", "pval": "pval_train"})
    te = test_fit.rename(columns={"rho": "rho_test", "pval": "pval_test"})
    # ANALYSIS_OK[join]: feature_index is unique within each of tr/te; constant-feature NaNs dropped after the join, not silently coerced.
    merged = tr.merge(te, on="feature_index", validate="one_to_one").dropna()

    # ANALYSIS_OK[positional-access]: 2x2 correlation matrix from .corr() -- [0,1] is always the train/test cross-term for a 2-column input.
    spearman_calib = merged[["rho_train", "rho_test"]].corr(method="spearman").iloc[0, 1]
    print(f"\nglobal calibration (Spearman rho of rho_train vs rho_test, {len(merged)} features): {spearman_calib:.4f}")

    top_train = merged.sort_values("pval_train").head(TOP_N).copy()
    top_train["same_sign"] = np.sign(top_train["rho_train"]) == np.sign(top_train["rho_test"])
    top_train["replicated"] = top_train["same_sign"] & (top_train["pval_test"] < 0.05)
    n_replicated = int(top_train["replicated"].sum())
    n_tested = len(top_train)
    print(f"train top-{TOP_N} hit replication in held-out test: {n_replicated}/{n_tested} "
          "(same sign AND test p<0.05)")
    print(top_train[["feature_index", "rho_train", "pval_train", "rho_test", "pval_test", "replicated"]].to_string(index=False))

    # ANALYSIS_OK[file-selection]: os.path.join(out_dir, ...) writes a named per-track output file, not a glob/latest read.
    merged.to_csv(os.path.join(out_dir, "holdout_calibration_all_features.csv.gz"), index=False, compression="gzip")
    # ANALYSIS_OK[file-selection]: os.path.join(out_dir, ...) writes a named per-track output file, not a glob/latest read.
    top_train.to_csv(os.path.join(out_dir, "holdout_hit_replication.csv"), index=False)

    track_summary = {
        "n_species_splittable": int(len(splittable_species)),
        "n_species_total": int(len(species_counts)),
        "n_train_rows": n_train_rows,
        "n_test_rows": n_test_rows,
        "spearman_calibration_train_vs_test": float(spearman_calib),
        "top_hit_replication_rate": n_replicated / n_tested if n_tested else None,
        "n_top_hits_tested": n_tested,
        "n_top_hits_replicated": n_replicated,
    }
    overall_summary[track] = track_summary
    # ANALYSIS_OK[file-selection]: provenance string names a fixed, just-written per-track output file, not a glob/latest read.
    register_value(f"{track}_holdout_replication_rate", track_summary["top_hit_replication_rate"], provenance=f"outputs/{track}/holdout_hit_replication.csv")
    # ANALYSIS_OK[file-selection]: provenance string names a fixed, just-written per-track output file, not a glob/latest read.
    register_value(f"{track}_holdout_calibration_spearman", track_summary["spearman_calibration_train_vs_test"], provenance=f"outputs/{track}/holdout_calibration_all_features.csv.gz")

with open(os.path.join(OUT_DIR, "holdout_summary.json"), "w") as f:
    json.dump(overall_summary, f, indent=2)
print("\n" + json.dumps(overall_summary, indent=2))
