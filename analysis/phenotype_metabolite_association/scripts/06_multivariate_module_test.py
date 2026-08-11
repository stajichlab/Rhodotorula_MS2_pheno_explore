#!/usr/bin/env python3
"""Multivariate/module-level test: is there a joint signal across many metabolite
features that predicts color phenotype, even though no single feature clears
02_corrected_correlation.py's univariate FDR test?

Method: PLS regression (X = rank-residualized 7,341 features, Y = rank-residualized
[L*, a*, b*], both adjusted for Species + Library Plate + sample_type exactly as in
02_corrected_correlation.py) evaluated by out-of-fold R^2 under GroupKFold
cross-validation, grouped by strain_id so a strain's C_*/SUP_* rows never split
across train/test (same non-independence concern as everywhere else in this
analysis). The number of PLS components is chosen by the CV itself (best mean
out-of-fold R^2 over a small grid), and -- critically -- the *entire* grid-search
procedure is repeated inside a strain-block permutation null (same design as
03_permutation_null.py) so the reported p-value already accounts for the
"pick the best of several component counts" selection, not just a single fixed
model. This directly answers whether a joint signal across features exists at all,
independent of which single feature it's carried by.
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

PHENOTYPE_COLS = ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
N_COMPONENTS_GRID = [2, 5, 10]
N_SPLITS = 4
N_PERM = int(os.environ.get("N_PERM", "300"))
RNG_SEED = 20260811

print("=" * 80)
print("06: MULTIVARIATE / MODULE-LEVEL TEST (PLS regression, CV R^2, permutation null)")
print("=" * 80)

design = pd.read_csv(os.path.join(OUT_DIR, "sample_design.csv"))
features = pd.read_csv(os.path.join(OUT_DIR, "features_cleaned.csv.gz"))
n = len(design)
rng = np.random.default_rng(RNG_SEED)

# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against
# Species/Library Plate/sample_type, identical to 02_corrected_correlation.py --
# see that script's docstring for the two bugs this fixes in the legacy pipeline.
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
# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type (or Plate/sample_type only within the single-species subset); same method as 02_corrected_correlation.py, documented in the module docstring.
X = residualize(franks)
X = (X - X.mean(axis=0)) / X.std(axis=0)  # PLS is scale-sensitive; features already rank-based

Y_ranks = np.column_stack([rankdata(design[c].values.astype(float)) for c in PHENOTYPE_COLS])
# ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type (or Plate/sample_type only within the single-species subset); same method as 02_corrected_correlation.py, documented in the module docstring.
Y = residualize(Y_ranks)
Y = (Y - Y.mean(axis=0)) / Y.std(axis=0)

strain_id = design["strain_id"].values

# --- strain-block permutation machinery, same design as 03_permutation_null.py ---
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


def best_cv_r2(Y_this):
    """Grid-search n_components by GroupKFold CV R^2; return the best mean R^2
    and the component count that achieved it. This whole procedure -- not a
    single fixed model -- is what gets permutation-tested below."""
    gkf = GroupKFold(n_splits=N_SPLITS)
    best = (-np.inf, None)
    for k in N_COMPONENTS_GRID:
        fold_r2 = []
        for train_idx, test_idx in gkf.split(X, Y_this, groups=strain_id):
            model = PLSRegression(n_components=k)
            model.fit(X[train_idx], Y_this[train_idx])
            pred = model.predict(X[test_idx])
            fold_r2.append(r2_score(Y_this[test_idx], pred, multioutput="uniform_average"))
        mean_r2 = float(np.mean(fold_r2))
        if mean_r2 > best[0]:
            best = (mean_r2, k)
    return best


print(f"\n[1/2] Observed CV R^2 (grid={N_COMPONENTS_GRID}, {N_SPLITS}-fold GroupKFold by strain)...")
t0 = time.time()
obs_r2, obs_k = best_cv_r2(Y)
print(f"  best mean CV R^2 = {obs_r2:.4f} at n_components={obs_k}  ({time.time()-t0:.1f}s)")

# per-phenotype breakdown at the winning component count, for interpretability
gkf = GroupKFold(n_splits=N_SPLITS)
per_pheno_r2 = {c: [] for c in PHENOTYPE_COLS}
for train_idx, test_idx in gkf.split(X, Y, groups=strain_id):
    model = PLSRegression(n_components=obs_k)
    model.fit(X[train_idx], Y[train_idx])
    pred = model.predict(X[test_idx])
    for j, c in enumerate(PHENOTYPE_COLS):
        per_pheno_r2[c].append(r2_score(Y[test_idx, j], pred[:, j]))
per_pheno_r2_mean = {c: float(np.mean(v)) for c, v in per_pheno_r2.items()}
print("  per-phenotype CV R^2 at winning n_components:")
for c, v in per_pheno_r2_mean.items():
    print(f"    {c:28} {v:7.4f}")

print(f"\n[2/2] Permutation null ({N_PERM} strain-block permutations, full grid-search repeated each time)...")
t0 = time.time()
null_r2 = np.empty(N_PERM)
for p in range(N_PERM):
    perm = permuted_row_order()
    Y_perm = Y[perm]
    null_r2[p], _ = best_cv_r2(Y_perm)
    if (p + 1) % max(1, N_PERM // 10) == 0:
        elapsed = time.time() - t0
        print(f"  {p+1}/{N_PERM} permutations done ({elapsed:.0f}s elapsed, "
              f"~{elapsed / (p + 1) * (N_PERM - p - 1):.0f}s remaining)")

perm_pval = (np.sum(null_r2 >= obs_r2) + 1) / (N_PERM + 1)
null_2p5, null_97p5 = np.percentile(null_r2, [2.5, 97.5])
print(f"\nobserved best CV R^2 = {obs_r2:.4f} (n_components={obs_k})")
print(f"null CV R^2: mean={null_r2.mean():.4f}  95% range=[{null_2p5:.4f}, {null_97p5:.4f}]")
print(f"permutation p-value = {perm_pval:.4f}")

pd.DataFrame({"null_r2": null_r2}).to_csv(os.path.join(OUT_DIR, "multivariate_permutation_null.csv"), index=False)

summary = {
    "n_components_grid": N_COMPONENTS_GRID,
    "n_splits": N_SPLITS,
    "n_permutations": N_PERM,
    "observed_best_cv_r2": obs_r2,
    "observed_best_n_components": obs_k,
    "per_phenotype_cv_r2_at_winning_k": per_pheno_r2_mean,
    "null_r2_mean": float(null_r2.mean()),
    "null_r2_2p5": float(null_2p5),
    "null_r2_97p5": float(null_97p5),
    "permutation_pvalue": float(perm_pval),
    "rng_seed": RNG_SEED,
}
with open(os.path.join(OUT_DIR, "multivariate_module_test_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n" + json.dumps(summary, indent=2))

register_value("multivariate_observed_cv_r2", obs_r2, provenance="outputs/multivariate_module_test_summary.json")
register_value("multivariate_permutation_pvalue", float(perm_pval), provenance="outputs/multivariate_module_test_summary.json")
