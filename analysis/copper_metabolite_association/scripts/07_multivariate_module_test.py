#!/usr/bin/env python3
"""Multivariate/module-level test, per track: is there a joint signal across many
metabolite features that predicts copper AUC, even though no single feature clears
03_corrected_correlation.py's univariate FDR test?

Method: PLS regression (X = rank-residualized features, y = rank-residualized
mean_auc_rate, both adjusted for the track's covariates exactly as in
03_corrected_correlation.py) evaluated by out-of-fold R^2 under GroupKFold CV,
grouped by strain_id (each track already has one row per strain -- no C/SUP pairing
within a track -- so this just prevents the same strain appearing in both grid-search
folds if a future version adds replicate rows; harmless no-op otherwise). Number of
PLS components chosen by CV itself over a small grid, and the *entire* grid-search is
repeated inside a permutation null (simple row permutation -- see 05_permutation_null.py
for why no block structure is needed here) so the reported p-value accounts for the
component-count selection, per the nested-design fix validated in
analysis/phenotype_metabolite_association/scripts/06_multivariate_module_test.py.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import r2_score
from sklearn.model_selection import GroupKFold

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
N_COMPONENTS_GRID = [2, 5, 10]
N_SPLITS = 4
N_PERM = int(os.environ.get("N_PERM", "300"))
RNG_SEED = 20260811

print("=" * 80)
print("07: MULTIVARIATE / MODULE-LEVEL TEST (PLS regression, CV R^2, permutation null)")
print("=" * 80)


def run_track(track, covariate_cols):
    print("\n" + "-" * 80)
    print(f"TRACK: {track}")
    print("-" * 80)
    out_dir = os.path.join(OUT_DIR, track)
    rng = np.random.default_rng(RNG_SEED)

    design = pd.read_csv(os.path.join(out_dir, "sample_design.csv"))
    features = pd.read_csv(os.path.join(out_dir, "features_cleaned.csv.gz"))
    n = len(design)

    Z_parts = [np.ones((n, 1))]
    for col in covariate_cols:
        dummies = pd.get_dummies(design[col].astype(str), drop_first=True)
        Z_parts.append(dummies.values.astype(float))
    Z = np.hstack(Z_parts)
    Z_pinv = np.linalg.pinv(Z)

    def residualize(M):
        return M - Z @ (Z_pinv @ M)

    fvals = features.values.astype(float)
    franks = np.apply_along_axis(rankdata, 0, fvals)
    # ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py.
    X = residualize(franks)
    X_std = X.std(axis=0)
    X_std[X_std < 1e-10] = np.nan  # constant features -> NaN column, dropped next
    X = (X - X.mean(axis=0)) / X_std
    keep_cols = ~np.isnan(X).any(axis=0)
    n_dropped = int((~keep_cols).sum())
    X = X[:, keep_cols]
    print(f"features: {fvals.shape[1]} total, {n_dropped} dropped (zero residual variance), {X.shape[1]} used")

    y_rank = rankdata(design["mean_auc_rate"].values.astype(float))
    # ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py.
    y = residualize(y_rank.reshape(-1, 1)).ravel()
    y = (y - y.mean()) / y.std()

    strain_id = design["strain_id"].values

    def best_cv_r2(y_this):
        gkf = GroupKFold(n_splits=N_SPLITS)
        best = (-np.inf, None)
        for k in N_COMPONENTS_GRID:
            fold_r2 = []
            for train_idx, test_idx in gkf.split(X, y_this, groups=strain_id):
                model = PLSRegression(n_components=k)
                model.fit(X[train_idx], y_this[train_idx])
                pred = model.predict(X[test_idx]).ravel()
                fold_r2.append(r2_score(y_this[test_idx], pred))
            mean_r2 = float(np.mean(fold_r2))
            if mean_r2 > best[0]:
                best = (mean_r2, k)
        return best

    print(f"\n[1/2] Observed CV R^2 (grid={N_COMPONENTS_GRID}, {N_SPLITS}-fold GroupKFold by strain)...")
    t0 = time.time()
    obs_r2, obs_k = best_cv_r2(y)
    print(f"  best mean CV R^2 = {obs_r2:.4f} at n_components={obs_k}  ({time.time()-t0:.1f}s)")

    print(f"\n[2/2] Permutation null ({N_PERM} row permutations, full grid-search repeated each time)...")
    t0 = time.time()
    null_r2 = np.empty(N_PERM)
    for p in range(N_PERM):
        perm = rng.permutation(n)
        null_r2[p], _ = best_cv_r2(y[perm])
        if (p + 1) % max(1, N_PERM // 5) == 0:
            elapsed = time.time() - t0
            print(f"  {p+1}/{N_PERM} ({elapsed:.0f}s elapsed, ~{elapsed / (p + 1) * (N_PERM - p - 1):.0f}s remaining)")

    perm_pval = (np.sum(null_r2 >= obs_r2) + 1) / (N_PERM + 1)
    null_2p5, null_97p5 = np.percentile(null_r2, [2.5, 97.5])
    print(f"\nobserved best CV R^2 = {obs_r2:.4f} (n_components={obs_k})")
    print(f"null CV R^2: mean={null_r2.mean():.4f}  95% range=[{null_2p5:.4f}, {null_97p5:.4f}]")
    print(f"permutation p-value = {perm_pval:.4f}")

    # ANALYSIS_OK[file-selection]: named per-track output file, not a glob/latest read.
    pd.DataFrame({"null_r2": null_r2}).to_csv(os.path.join(out_dir, "multivariate_permutation_null.csv"), index=False)

    return {
        "n_samples": int(n),
        "n_features_used": int(X.shape[1]),
        "n_components_grid": N_COMPONENTS_GRID,
        "n_splits": N_SPLITS,
        "n_permutations": N_PERM,
        "observed_best_cv_r2": float(obs_r2),
        "observed_best_n_components": obs_k,
        "null_r2_mean": float(null_r2.mean()),
        "null_r2_2p5": float(null_2p5),
        "null_r2_97p5": float(null_97p5),
        "permutation_pvalue": float(perm_pval),
        "rng_seed": RNG_SEED,
    }


overall = {}
for track, covariate_cols in TRACK_COVARIATES.items():
    overall[track] = run_track(track, covariate_cols)
    register_value(f"{track}_multivariate_cv_r2", overall[track]["observed_best_cv_r2"], provenance=f"outputs/{track}/multivariate_module_test_summary.json")
    register_value(f"{track}_multivariate_perm_pvalue", overall[track]["permutation_pvalue"], provenance=f"outputs/{track}/multivariate_module_test_summary.json")

with open(os.path.join(OUT_DIR, "multivariate_module_test_summary.json"), "w") as f:
    json.dump(overall, f, indent=2)
print("\n" + json.dumps(overall, indent=2))
