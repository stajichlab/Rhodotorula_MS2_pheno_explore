#!/usr/bin/env python3
"""Cross-reference this analysis's cell-track multivariate metabolite signal
(F-004/F-005: PLS+RF both permutation-significant for the cell fraction, see
COPPER_METABOLITE_ASSOCIATION.md) against the independent genome-side copper signal
found in the sibling Rhodotorula_Rodeo repo (species-corrected whole-proteome
methionine usage correlates with copper AUC, `.living/findings/metal-resistance.md`
F-002 there). This is the actual "connect metabolites to secreted-protein/genome
predictions" question the user originally asked about (2026-08-11) -- everything
before this script was building the two independent signals; this is where they meet.

Question: do strains whose CELL-fraction metabolite profile better predicts high
copper AUC (out-of-fold PLS/RF prediction, in species-corrected rank-residual space)
also carry higher species-corrected genome methionine usage?

Method:
1. Refit the cell track's winning PLS (n_components=10) and Random Forest
   (max_depth=None) models under the same GroupKFold CV as
   07/08_*.py, but this time save the OUT-OF-FOLD PREDICTED y_resid per strain
   (an honest, non-overfit "how much does this strain's metabolite profile look
   like a high-copper-AUC profile" score) instead of just the aggregate R^2.
2. Join to Rhodotorula_Rodeo's met_vs_copper_auc_data.csv by normalized strain name
   (SAMPLE_NAME <-> LOCUSTAG; the cell track's own `strain_id` column uses
   ATTRIBUTE_ID_1, which drops the TFCN_/NRRL_ prefix that Rodeo's LOCUSTAG keeps --
   SAMPLE_NAME is the correct join key here, confirmed against Rodeo's own
   01b_build_continuous_cohort.py join logic).
3. Spearman-correlate the out-of-fold metabolite prediction against Rodeo's
   species-corrected Met-fraction residual, across the common strain set, with a
   permutation null.

This does NOT re-derive or duplicate Rodeo's data -- it reads
`metal_resistance/results/figures/met_vs_copper_auc_data.csv` in place from that
repo, the same read-only cross-repo pattern already used when building that file's
underlying cohort join.
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "outputs")
RODEO_MET_CSV = os.path.join(
    ROOT, "..", "..", "Rhodotorula_Rodeo", "metal_resistance", "results", "figures",
    "met_vs_copper_auc_data.csv",
)

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

TRACK = "cell"  # only the cell track had cross-model-class-corroborated signal (F-005)
COVARIATE_COLS = ["species", "Library Plate"]
PLS_N_COMPONENTS = 10  # winning k from 07_multivariate_module_test.py's cell-track grid search
RF_MAX_DEPTH = None  # winning depth from 08_random_forest_module_test.py's cell-track grid search
RF_N_ESTIMATORS = 150
N_SPLITS = 4
N_PERM = 5000
RNG_SEED = 20260811


def norm(s):
    return re.sub(r"[^A-Za-z0-9]", "", str(s)).upper()


print("=" * 80)
print("09: CROSS-REFERENCE cell-track metabolite signal vs. Rhodotorula_Rodeo genome Met signal")
print("=" * 80)

if not os.path.exists(RODEO_MET_CSV):
    raise FileNotFoundError(
        f"{RODEO_MET_CSV} not found -- run Rhodotorula_Rodeo's "
        "metal_resistance/scripts/04c_plot_met_correlation.py first (see that repo's "
        "PLAN_metal_resistance.md, Phase 1b)."
    )

out_dir = os.path.join(OUT_DIR, TRACK)
design = pd.read_csv(os.path.join(out_dir, "sample_design.csv"))
features = pd.read_csv(os.path.join(out_dir, "features_cleaned.csv.gz"))
n = len(design)
print(f"cell-track strains: {n}")

Z_parts = [np.ones((n, 1))]
for col in COVARIATE_COLS:
    dummies = pd.get_dummies(design[col].astype(str), drop_first=True)
    Z_parts.append(dummies.values.astype(float))
Z = np.hstack(Z_parts)
Z_pinv = np.linalg.pinv(Z)


def residualize(M):
    return M - Z @ (Z_pinv @ M)


fvals = features.values.astype(float)
franks = np.apply_along_axis(rankdata, 0, fvals)
# ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py / 07/08_*.py.
X_full = residualize(franks)
X_std = X_full.std(axis=0)
X_std[X_std < 1e-10] = np.nan
X_full = (X_full - X_full.mean(axis=0)) / X_std
keep_cols = ~np.isnan(X_full).any(axis=0)
X = X_full[:, keep_cols]
print(f"features used: {X.shape[1]} (of {fvals.shape[1]})")

y_rank = rankdata(design["mean_auc_rate"].values.astype(float))
# ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py / 07/08_*.py.
y = residualize(y_rank.reshape(-1, 1)).ravel()
y = (y - y.mean()) / y.std()

strain_id = design["strain_id"].values

print(f"\n[1/3] Out-of-fold predictions: PLS (n_components={PLS_N_COMPONENTS}), "
      f"RF (max_depth={RF_MAX_DEPTH}), {N_SPLITS}-fold GroupKFold by strain...")
gkf = GroupKFold(n_splits=N_SPLITS)
pls_pred = np.full(n, np.nan)
rf_pred = np.full(n, np.nan)
for train_idx, test_idx in gkf.split(X, y, groups=strain_id):
    pls = PLSRegression(n_components=PLS_N_COMPONENTS)
    pls.fit(X[train_idx], y[train_idx])
    pls_pred[test_idx] = pls.predict(X[test_idx]).ravel()

    rf = RandomForestRegressor(
        n_estimators=RF_N_ESTIMATORS, max_depth=RF_MAX_DEPTH, max_features="sqrt",
        n_jobs=-1, random_state=RNG_SEED,
    )
    rf.fit(X[train_idx], y[train_idx])
    rf_pred[test_idx] = rf.predict(X[test_idx])

# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert not np.isnan(pls_pred).any() and not np.isnan(rf_pred).any(), "every row must get exactly one out-of-fold prediction from GroupKFold"

pred_df = design[["SAMPLE_NAME", "strain_id", "species", "mean_auc_rate"]].copy()
pred_df["pls_oof_pred"] = pls_pred
pred_df["rf_oof_pred"] = rf_pred
pred_df["norm_name"] = pred_df["SAMPLE_NAME"].map(norm)

print("\n[2/3] Joining to Rhodotorula_Rodeo genome Met-usage data...")
met = pd.read_csv(RODEO_MET_CSV)
met["norm_name"] = met["LOCUSTAG"].map(norm)
# ANALYSIS_OK[join]: many_to_one is expected only if pred_df's norm_name is
# unique, which sample_design.csv construction (01_prepare_data.py) guarantees
# (strain_id/SAMPLE_NAME uniqueness asserted there); validated explicitly here too.
assert pred_df["norm_name"].is_unique, "cell-track SAMPLE_NAME must be unique after normalization"
assert met["norm_name"].is_unique, "Rodeo met_vs_copper_auc_data LOCUSTAG must be unique after normalization"

merged = pred_df.merge(
    met[["norm_name", "met_fraction", "met_fraction_rank_residual_species", "copper_auc_rank_residual_species"]],
    on="norm_name", how="inner", validate="one_to_one",
)
n_common = len(merged)
print(f"strains present in both analyses: {n_common} / {n} cell-track, {len(met)} Rodeo genome-Met")
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert n_common > 30, "too few common strains for a meaningful cross-reference"

# sanity check: mean_auc_rate should agree between the two independently-loaded copies
max_auc_diff = float((merged["mean_auc_rate"] - met.set_index("norm_name").loc[merged["norm_name"], "mean_auc_rate"].values).abs().max())
print(f"sanity check: max |mean_auc_rate difference| between independently-loaded copies = {max_auc_diff:.2e} (should be ~0)")
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert max_auc_diff < 1e-6, "mean_auc_rate mismatch between the two repos' copies of the same phenotype -- join is wrong"

print(f"\n[3/3] Correlating out-of-fold metabolite predictions against genome Met residual "
      f"(n={n_common} common strains, {N_PERM} permutations)...")


def perm_test(a, b, n_perm=N_PERM, seed=RNG_SEED):
    rng = np.random.default_rng(seed)
    obs_rho, obs_p_parametric = spearmanr(a, b)
    count = 0
    for _ in range(n_perm):
        b_perm = rng.permutation(b)
        rho, _ = spearmanr(a, b_perm)
        if abs(rho) >= abs(obs_rho):
            count += 1
    perm_p = (count + 1) / (n_perm + 1)
    return float(obs_rho), float(obs_p_parametric), float(perm_p)


results = {}
for model_col, model_name in [("pls_oof_pred", "PLS"), ("rf_oof_pred", "RandomForest")]:
    rho, p_param, p_perm = perm_test(merged[model_col].values, merged["met_fraction_rank_residual_species"].values)
    print(f"  {model_name} out-of-fold prediction vs. genome Met residual: "
          f"rho={rho:+.4f}, parametric p={p_param:.4g}, permutation p={p_perm:.4f}")
    results[model_name] = {"rho": rho, "parametric_p": p_param, "permutation_p": p_perm}

merged.to_csv(os.path.join(out_dir, "cross_reference_genome_met.csv"), index=False)
summary = {
    "track": TRACK,
    "n_cell_track_strains": int(n),
    "n_rodeo_met_strains": int(len(met)),
    "n_common_strains": int(n_common),
    "n_permutations": N_PERM,
    "rng_seed": RNG_SEED,
    "correlations": results,
}
with open(os.path.join(OUT_DIR, "cross_reference_genome_met_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n✓ wrote {out_dir}/cross_reference_genome_met.csv, "
      f"{OUT_DIR}/cross_reference_genome_met_summary.json")
print("\n" + json.dumps(summary, indent=2))

for model_name, r in results.items():
    register_value(f"cross_ref_{model_name}_rho", r["rho"], provenance="outputs/cross_reference_genome_met_summary.json")
    register_value(f"cross_ref_{model_name}_permutation_p", r["permutation_p"], provenance="outputs/cross_reference_genome_met_summary.json")
