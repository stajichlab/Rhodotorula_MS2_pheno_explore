#!/usr/bin/env python3
"""Non-linear multivariate/module-level test, per track: does a Random Forest find
a joint signal that 07_multivariate_module_test.py's (linear) PLS missed?

Same design as script 07 for direct comparability, adapted from
analysis/phenotype_metabolite_association/scripts/08_random_forest_module_test.py:
rank-residualized features (X) and rank-residualized copper AUC (y), adjusted for
the track's covariates; out-of-fold R^2 under GroupKFold (by strain_id); a small
max_depth grid chosen by CV; the entire grid-search re-run inside a permutation null
(simple row permutation, no block structure needed -- see 05_permutation_null.py).
`max_features='sqrt'` for tractable runtime at ~7,300 features, same rationale as
the color-phenotype analysis.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.ensemble import RandomForestRegressor
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
MAX_DEPTH_GRID = [6, None]
N_ESTIMATORS = 150
N_SPLITS = 4
N_PERM = int(os.environ.get("N_PERM", "200"))
RNG_SEED = 20260811

print("=" * 80)
print("08: RANDOM FOREST MODULE-LEVEL TEST (non-linear, CV R^2, permutation null)")
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
    X = residualize(franks).astype(np.float32)  # float32: RF fit cost, no precision loss that matters for ranks

    y_rank = rankdata(design["mean_auc_rate"].values.astype(float))
    # ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py.
    y = residualize(y_rank.reshape(-1, 1)).ravel()

    strain_id = design["strain_id"].values

    def best_cv_r2(y_this, seed):
        gkf = GroupKFold(n_splits=N_SPLITS)
        best = (-np.inf, None)
        for depth in MAX_DEPTH_GRID:
            fold_r2 = []
            for train_idx, test_idx in gkf.split(X, y_this, groups=strain_id):
                model = RandomForestRegressor(
                    n_estimators=N_ESTIMATORS, max_depth=depth, max_features="sqrt",
                    n_jobs=-1, random_state=seed,
                )
                model.fit(X[train_idx], y_this[train_idx])
                pred = model.predict(X[test_idx])
                fold_r2.append(r2_score(y_this[test_idx], pred))
            mean_r2 = float(np.mean(fold_r2))
            if mean_r2 > best[0]:
                best = (mean_r2, depth)
        return best

    print(f"\n[1/3] Observed CV R^2 (grid max_depth={MAX_DEPTH_GRID}, n_estimators={N_ESTIMATORS}, "
          f"{N_SPLITS}-fold GroupKFold by strain)...")
    t0 = time.time()
    obs_r2, obs_depth = best_cv_r2(y, seed=RNG_SEED)
    print(f"  best mean CV R^2 = {obs_r2:.4f} at max_depth={obs_depth}  ({time.time()-t0:.1f}s)")

    print("\n[2/3] Feature importances at the winning depth...")
    gkf = GroupKFold(n_splits=N_SPLITS)
    importances = np.zeros(X.shape[1])
    for train_idx, test_idx in gkf.split(X, y, groups=strain_id):
        model = RandomForestRegressor(
            n_estimators=N_ESTIMATORS, max_depth=obs_depth, max_features="sqrt",
            n_jobs=-1, random_state=RNG_SEED,
        )
        model.fit(X[train_idx], y[train_idx])
        importances += model.feature_importances_
    importances /= N_SPLITS

    top_importance = pd.DataFrame({"feature_index": np.arange(X.shape[1]), "mean_importance": importances})
    top_importance = top_importance.sort_values("mean_importance", ascending=False)
    # ANALYSIS_OK[file-selection]: named per-track output file, not a glob/latest read.
    top_importance.to_csv(os.path.join(out_dir, "random_forest_feature_importances.csv"), index=False)
    print("  top 10 features by mean CV feature_importances_ (descriptive only -- see permutation "
          "p-value below before treating these as candidates):")
    print(top_importance.head(10).to_string(index=False))

    print(f"\n[3/3] Permutation null ({N_PERM} row permutations, full grid-search repeated each time)...")
    t0 = time.time()
    null_r2 = np.empty(N_PERM)
    for p in range(N_PERM):
        perm = rng.permutation(n)
        null_r2[p], _ = best_cv_r2(y[perm], seed=RNG_SEED + p + 1)
        if (p + 1) % max(1, N_PERM // 5) == 0:
            elapsed = time.time() - t0
            print(f"  {p+1}/{N_PERM} ({elapsed:.0f}s elapsed, ~{elapsed / (p + 1) * (N_PERM - p - 1):.0f}s remaining)")

    perm_pval = (np.sum(null_r2 >= obs_r2) + 1) / (N_PERM + 1)
    null_2p5, null_97p5 = np.percentile(null_r2, [2.5, 97.5])
    print(f"\nobserved best CV R^2 = {obs_r2:.4f} (max_depth={obs_depth})")
    print(f"null CV R^2: mean={null_r2.mean():.4f}  95% range=[{null_2p5:.4f}, {null_97p5:.4f}]")
    print(f"permutation p-value = {perm_pval:.4f}")

    # ANALYSIS_OK[file-selection]: named per-track output file, not a glob/latest read.
    pd.DataFrame({"null_r2": null_r2}).to_csv(os.path.join(out_dir, "random_forest_permutation_null.csv"), index=False)

    return {
        "n_samples": int(n),
        "max_depth_grid": [str(d) for d in MAX_DEPTH_GRID],
        "n_estimators": N_ESTIMATORS,
        "max_features": "sqrt",
        "n_splits": N_SPLITS,
        "n_permutations": N_PERM,
        "observed_best_cv_r2": float(obs_r2),
        "observed_best_max_depth": str(obs_depth),
        "null_r2_mean": float(null_r2.mean()),
        "null_r2_2p5": float(null_2p5),
        "null_r2_97p5": float(null_97p5),
        "permutation_pvalue": float(perm_pval),
        "rng_seed": RNG_SEED,
    }


overall = {}
for track, covariate_cols in TRACK_COVARIATES.items():
    overall[track] = run_track(track, covariate_cols)
    register_value(f"{track}_rf_cv_r2", overall[track]["observed_best_cv_r2"], provenance=f"outputs/{track}/random_forest_module_test_summary.json")
    register_value(f"{track}_rf_perm_pvalue", overall[track]["permutation_pvalue"], provenance=f"outputs/{track}/random_forest_module_test_summary.json")

with open(os.path.join(OUT_DIR, "random_forest_module_test_summary.json"), "w") as f:
    json.dump(overall, f, indent=2)
print("\n" + json.dumps(overall, indent=2))
