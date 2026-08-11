#!/usr/bin/env python3
"""Non-linear multivariate/module-level test: does a Random Forest, which can
capture feature interactions and non-linear metabolite-phenotype relationships
that PLS (linear) cannot, find a joint signal that 06_multivariate_module_test.py
missed?

Same design as script 06 for direct comparability: rank-residualized features (X)
and rank-residualized phenotype (Y = [L*, a*, b*]), both adjusted for Species +
Library Plate + sample_type; out-of-fold R^2 under GroupKFold (grouped by
strain_id); a small hyperparameter grid (max_depth) chosen by the CV itself; and
the *entire* grid-search re-run inside a strain-block permutation null so the
p-value accounts for the max_depth selection, not just one fixed model.
`max_features='sqrt'` is required for tractable runtime at 7,341 features (the
sklearn default of considering all features at every split is >100x slower here
and buys no accuracy at this n/p ratio) -- confirmed empirically before writing
this script.
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

PHENOTYPE_COLS = ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
MAX_DEPTH_GRID = [6, None]
N_ESTIMATORS = 150
N_SPLITS = 4
N_PERM = int(os.environ.get("N_PERM", "200"))
RNG_SEED = 20260811

print("=" * 80)
print("08: RANDOM FOREST MODULE-LEVEL TEST (non-linear, CV R^2, permutation null)")
print("=" * 80)

design = pd.read_csv(os.path.join(OUT_DIR, "sample_design.csv"))
features = pd.read_csv(os.path.join(OUT_DIR, "features_cleaned.csv.gz"))
n = len(design)
rng = np.random.default_rng(RNG_SEED)

# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against
# Species/Library Plate/sample_type, identical to 02_corrected_correlation.py and
# 06_multivariate_module_test.py -- see those docstrings for the two bugs this
# fixes relative to the legacy pipeline.
Z_parts = [np.ones((n, 1))]
for col in ("species", "Library Plate", "sample_type"):
    dummies = pd.get_dummies(design[col].astype(str), drop_first=True)
    Z_parts.append(dummies.values.astype(float))
Z = np.hstack(Z_parts)
Z_pinv = np.linalg.pinv(Z)


def residualize(M):
    return M - Z @ (Z_pinv @ M)


fvals = features.values.astype(float)
franks = np.apply_along_axis(rankdata, 0, fvals)
# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type, identical to 02_corrected_correlation.py and 06_multivariate_module_test.py; documented in the module docstring.
X = residualize(franks).astype(np.float32)  # float32: RF fit cost, no precision loss that matters for ranks

Y_ranks = np.column_stack([rankdata(design[c].values.astype(float)) for c in PHENOTYPE_COLS])
# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type, identical to 02_corrected_correlation.py and 06_multivariate_module_test.py; documented in the module docstring.
Y = residualize(Y_ranks)

strain_id = design["strain_id"].values

# --- strain-block permutation machinery, identical design to 03/06 ---
type_of = design["sample_type"].values
strain_rows = {}
for i, s in enumerate(strain_id):
    strain_rows.setdefault(s, {})[type_of[i]] = i
paired = [s for s, d in strain_rows.items() if set(d) == {"C", "SUP"}]
c_only = [s for s, d in strain_rows.items() if set(d) == {"C"}]
sup_only = [s for s, d in strain_rows.items() if set(d) == {"SUP"}]


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


def best_cv_r2(Y_this, seed):
    """Grid-search max_depth by GroupKFold CV R^2; return the best mean R^2 and
    the depth that achieved it -- this whole procedure is what gets
    permutation-tested, not a single fixed model (same rationale as script 06)."""
    gkf = GroupKFold(n_splits=N_SPLITS)
    best = (-np.inf, None)
    for depth in MAX_DEPTH_GRID:
        fold_r2 = []
        for train_idx, test_idx in gkf.split(X, Y_this, groups=strain_id):
            model = RandomForestRegressor(
                n_estimators=N_ESTIMATORS, max_depth=depth, max_features="sqrt",
                n_jobs=-1, random_state=seed,
            )
            model.fit(X[train_idx], Y_this[train_idx])
            pred = model.predict(X[test_idx])
            fold_r2.append(r2_score(Y_this[test_idx], pred, multioutput="uniform_average"))
        mean_r2 = float(np.mean(fold_r2))
        if mean_r2 > best[0]:
            best = (mean_r2, depth)
    return best


print(f"\n[1/3] Observed CV R^2 (grid max_depth={MAX_DEPTH_GRID}, n_estimators={N_ESTIMATORS}, "
      f"{N_SPLITS}-fold GroupKFold by strain)...")
t0 = time.time()
obs_r2, obs_depth = best_cv_r2(Y, seed=RNG_SEED)
print(f"  best mean CV R^2 = {obs_r2:.4f} at max_depth={obs_depth}  ({time.time()-t0:.1f}s)")

print("\n[2/3] Per-phenotype breakdown and feature importances at the winning depth...")
gkf = GroupKFold(n_splits=N_SPLITS)
per_pheno_r2 = {c: [] for c in PHENOTYPE_COLS}
importances = np.zeros(X.shape[1])
for train_idx, test_idx in gkf.split(X, Y, groups=strain_id):
    model = RandomForestRegressor(
        n_estimators=N_ESTIMATORS, max_depth=obs_depth, max_features="sqrt",
        n_jobs=-1, random_state=RNG_SEED,
    )
    model.fit(X[train_idx], Y[train_idx])
    pred = model.predict(X[test_idx])
    importances += model.feature_importances_
    for j, c in enumerate(PHENOTYPE_COLS):
        per_pheno_r2[c].append(r2_score(Y[test_idx, j], pred[:, j]))
per_pheno_r2_mean = {c: float(np.mean(v)) for c, v in per_pheno_r2.items()}
importances /= N_SPLITS
print("  per-phenotype CV R^2 at winning max_depth:")
for c, v in per_pheno_r2_mean.items():
    print(f"    {c:28} {v:7.4f}")

top_importance = pd.DataFrame({"feature_index": np.arange(X.shape[1]), "mean_importance": importances})
top_importance = top_importance.sort_values("mean_importance", ascending=False)
top_importance.to_csv(os.path.join(OUT_DIR, "random_forest_feature_importances.csv"), index=False)
print("\n  top 10 features by mean CV feature_importances_ (descriptive only -- see permutation p-value below "
      "before treating these as candidates):")
print(top_importance.head(10).to_string(index=False))

print(f"\n[3/3] Permutation null ({N_PERM} strain-block permutations, full grid-search repeated each time)...")
t0 = time.time()
null_r2 = np.empty(N_PERM)
for p in range(N_PERM):
    perm = permuted_row_order()
    Y_perm = Y[perm]
    null_r2[p], _ = best_cv_r2(Y_perm, seed=RNG_SEED + p + 1)
    if (p + 1) % max(1, N_PERM // 10) == 0:
        elapsed = time.time() - t0
        print(f"  {p+1}/{N_PERM} permutations done ({elapsed:.0f}s elapsed, "
              f"~{elapsed / (p + 1) * (N_PERM - p - 1):.0f}s remaining)")

perm_pval = (np.sum(null_r2 >= obs_r2) + 1) / (N_PERM + 1)
null_2p5, null_97p5 = np.percentile(null_r2, [2.5, 97.5])
print(f"\nobserved best CV R^2 = {obs_r2:.4f} (max_depth={obs_depth})")
print(f"null CV R^2: mean={null_r2.mean():.4f}  95% range=[{null_2p5:.4f}, {null_97p5:.4f}]")
print(f"permutation p-value = {perm_pval:.4f}")

pd.DataFrame({"null_r2": null_r2}).to_csv(os.path.join(OUT_DIR, "random_forest_permutation_null.csv"), index=False)

summary = {
    "max_depth_grid": [str(d) for d in MAX_DEPTH_GRID],
    "n_estimators": N_ESTIMATORS,
    "max_features": "sqrt",
    "n_splits": N_SPLITS,
    "n_permutations": N_PERM,
    "observed_best_cv_r2": obs_r2,
    "observed_best_max_depth": str(obs_depth),
    "per_phenotype_cv_r2_at_winning_depth": per_pheno_r2_mean,
    "null_r2_mean": float(null_r2.mean()),
    "null_r2_2p5": float(null_2p5),
    "null_r2_97p5": float(null_97p5),
    "permutation_pvalue": float(perm_pval),
    "rng_seed": RNG_SEED,
}
with open(os.path.join(OUT_DIR, "random_forest_module_test_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n" + json.dumps(summary, indent=2))

register_value("rf_observed_cv_r2", obs_r2, provenance="outputs/random_forest_module_test_summary.json")
register_value("rf_permutation_pvalue", float(perm_pval), provenance="outputs/random_forest_module_test_summary.json")
