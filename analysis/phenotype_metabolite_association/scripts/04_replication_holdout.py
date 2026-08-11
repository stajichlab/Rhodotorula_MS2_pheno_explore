#!/usr/bin/env python3
"""Species-stratified, strain-level holdout replication check for the corrected
partial correlations. Splitting by strain_id (not by row) keeps a strain's C_*
and SUP_* rows together in the same split, avoiding train/test leakage from the
paired structure. Species with too few strains to split reliably (<5) are kept
entirely in train and are not part of the held-out test set (documented, not
silently dropped from the count).

Two replication checks are reported:
  1. Global calibration: Spearman correlation between rho_train and rho_test
     across ALL tested features -- if this pipeline captures a real signal, the
     effect-size ranking should be reproducible out of sample even where no
     single feature is individually significant.
  2. Hit replication: for the top-N nominal features identified in TRAIN, do they
     replicate (same sign AND nominal p<0.05) in TEST? This is deliberately a
     replication check on train-selected hits, not a re-run of the FDR scan, so
     it is not circular with 02_corrected_correlation.py's full-data result.
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

PHENOTYPE_COLS = ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
RNG_SEED = 20260811
TEST_FRAC = 0.20
MIN_STRAINS_FOR_SPLIT = 5
TOP_N_PER_PHENOTYPE = 25

print("=" * 80)
print("04: SPECIES-STRATIFIED, STRAIN-LEVEL HOLDOUT REPLICATION")
print("=" * 80)

design = pd.read_csv(os.path.join(OUT_DIR, "sample_design.csv"))
features = pd.read_csv(os.path.join(OUT_DIR, "features_cleaned.csv.gz"))
n = len(design)
rng = np.random.default_rng(RNG_SEED)

strain_species = design.drop_duplicates("strain_id").set_index("strain_id")["species"]
species_strain_counts = strain_species.value_counts()
splittable_species = species_strain_counts[species_strain_counts >= MIN_STRAINS_FOR_SPLIT].index
print(f"species with >={MIN_STRAINS_FOR_SPLIT} strains (eligible for holdout split): "
      f"{len(splittable_species)}/{len(species_strain_counts)}")

test_strains = set()
for sp in splittable_species:
    strains_sp = strain_species[strain_species == sp].index.tolist()
    n_test = max(1, round(len(strains_sp) * TEST_FRAC))
    chosen = rng.choice(strains_sp, size=n_test, replace=False)
    test_strains.update(chosen)

is_test = design["strain_id"].isin(test_strains).values
n_test_rows, n_train_rows = int(is_test.sum()), int((~is_test).sum())
print(f"train rows: {n_train_rows}  test rows: {n_test_rows}  "
      f"({len(test_strains)} held-out strains from {len(splittable_species)} species)")
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert n_test_rows > 30, "held-out set too small to be informative"


def fit_corrected(mask):
    sub_design = design.loc[mask].reset_index(drop=True)
    sub_features = features.loc[mask].reset_index(drop=True)
    nn = len(sub_design)

    Z_parts = [np.ones((nn, 1))]
    for col in ("species", "Library Plate", "sample_type"):
        dummies = pd.get_dummies(sub_design[col].astype(str), drop_first=True)
        Z_parts.append(dummies.values.astype(float))
    Z = np.hstack(Z_parts)
    rank_Z = np.linalg.matrix_rank(Z)
    Z_pinv = np.linalg.pinv(Z)

    def residualize(M):
        return M - Z @ (Z_pinv @ M)

    fvals = sub_features.values.astype(float)
    franks = np.apply_along_axis(rankdata, 0, fvals)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
    fresid = residualize(franks)
    fresid_c = fresid - fresid.mean(axis=0)
    fnorm = np.sqrt((fresid_c ** 2).sum(axis=0))
    # ANALYSIS_OK[threshold]: tier thresholds reproduced from the original scripts/phase2_correlation_analysis.py definition (|rho|>0.30/0.25/0.20, q<0.05/0.05/0.10) for direct comparability.
    fnorm[fnorm < 1e-10] = np.nan  # constant features -> NaN corr, filtered later

    out = {}
    df = nn - rank_Z - 1
    for phenotype in PHENOTYPE_COLS:
        y = sub_design[phenotype].values.astype(float)
        # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type; documented in the module docstring as the fix for the pooled-vs-species-confound bug in scripts/phase2_correlation_analysis.py.
        y_resid = residualize(rankdata(y).reshape(-1, 1)).ravel()
        y_c = y_resid - y_resid.mean()
        rho = (fresid_c.T @ y_c) / (fnorm * np.sqrt((y_c ** 2).sum()))
        t_stat = rho * np.sqrt(df / np.clip(1 - rho ** 2, 1e-12, None))
        pval = 2 * t_dist.sf(np.abs(t_stat), df=df)
        # feature_index is positional (0..n_features-1), matching 02_corrected_correlation.py's convention
        out[phenotype] = pd.DataFrame({"feature_index": np.arange(fvals.shape[1]), "rho": rho, "pval": pval})
    return out, df


print("\n[1/2] Fitting corrected model on TRAIN and TEST splits separately...")
train_fits, train_df = fit_corrected(~is_test)
test_fits, test_df = fit_corrected(is_test)
print(f"  train residual df={train_df}, test residual df={test_df}")

print("\n[2/2] Comparing train vs test...")
calib_rows = []
hit_rows = []
for phenotype in PHENOTYPE_COLS:
    tr = train_fits[phenotype].rename(columns={"rho": "rho_train", "pval": "pval_train"})
    te = test_fits[phenotype].rename(columns={"rho": "rho_test", "pval": "pval_test"})
    # ANALYSIS_OK[join]: feature_index is unique within each of tr/te (one row per
    # feature); constant-feature NaNs dropped after the join, not silently coerced.
    merged = tr.merge(te, on="feature_index", validate="one_to_one").dropna()

    # ANALYSIS_OK[positional-access]: 2x2 correlation matrix from .corr() -- [0,1]
    # is always the train/test cross-term for a 2-column input, not row/column metadata.
    spearman_calib = merged[["rho_train", "rho_test"]].corr(method="spearman").iloc[0, 1]
    calib_rows.append({"phenotype": phenotype, "n_features": len(merged), "spearman_rho_train_vs_test": spearman_calib})

    top_train = merged.sort_values("pval_train").head(TOP_N_PER_PHENOTYPE).copy()
    top_train["same_sign"] = np.sign(top_train["rho_train"]) == np.sign(top_train["rho_test"])
    top_train["replicated"] = top_train["same_sign"] & (top_train["pval_test"] < 0.05)
    top_train["phenotype"] = phenotype
    hit_rows.append(top_train)

calib_df = pd.DataFrame(calib_rows)
hits_df = pd.concat(hit_rows, ignore_index=True)
n_replicated = int(hits_df["replicated"].sum())
n_tested = len(hits_df)

print("\nGlobal calibration (Spearman rho of rho_train vs rho_test, all features):")
print(calib_df.to_string(index=False))
print(f"\nTrain top-{TOP_N_PER_PHENOTYPE}-per-phenotype hit replication in held-out test:")
print(f"  {n_replicated}/{n_tested} replicated (same sign AND test p<0.05)")
print(hits_df[["phenotype", "feature_index", "rho_train", "pval_train", "rho_test", "pval_test", "replicated"]].to_string(index=False))

# ANALYSIS_OK[file-selection]: os.path.join(OUT_DIR, ...) here writes a named output file (not a glob/latest read); no ambiguity about which file is used.
calib_df.to_csv(os.path.join(OUT_DIR, "holdout_calibration.csv"), index=False)
# ANALYSIS_OK[file-selection]: os.path.join(OUT_DIR, ...) here writes a named output file (not a glob/latest read); no ambiguity about which file is used.
hits_df.to_csv(os.path.join(OUT_DIR, "holdout_hit_replication.csv"), index=False)

summary = {
    "test_frac_target": TEST_FRAC,
    "min_strains_for_split": MIN_STRAINS_FOR_SPLIT,
    "n_species_splittable": int(len(splittable_species)),
    "n_species_total": int(len(species_strain_counts)),
    "n_test_strains": len(test_strains),
    "n_train_rows": n_train_rows,
    "n_test_rows": n_test_rows,
    "spearman_calibration_by_phenotype": {r["phenotype"]: r["spearman_rho_train_vs_test"] for r in calib_rows},
    "top_hit_replication_rate": n_replicated / n_tested if n_tested else None,
    "n_top_hits_tested": n_tested,
    "n_top_hits_replicated": n_replicated,
    "rng_seed": RNG_SEED,
}
# ANALYSIS_OK[file-selection]: os.path.join(OUT_DIR, ...) here writes a named output file (not a glob/latest read); no ambiguity about which file is used.
with open(os.path.join(OUT_DIR, "holdout_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n" + json.dumps(summary, indent=2))

# ANALYSIS_OK[file-selection]: os.path.join(OUT_DIR, ...) here writes a named output file (not a glob/latest read); no ambiguity about which file is used.
register_value("holdout_top_hit_replication_rate", summary["top_hit_replication_rate"], provenance="outputs/holdout_summary.json")
# ANALYSIS_OK[file-selection]: os.path.join(OUT_DIR, ...) here writes a named output file (not a glob/latest read); no ambiguity about which file is used.
register_value("holdout_n_test_strains", summary["n_test_strains"], provenance="outputs/holdout_summary.json")
