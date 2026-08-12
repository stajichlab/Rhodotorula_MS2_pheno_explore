#!/usr/bin/env python3
"""Power analysis, run BEFORE the permutation/holdout/multivariate steps -- explicit
lesson carried over from analysis/phenotype_metabolite_association/, where the power
analysis was written last and only then revealed the whole "five converging null
checks" result was uninformative below rho~0.34. Answering "could this design even
detect a real effect" up front changes how much weight the later null results should
carry.

Method: identical simulation approach to
analysis/phenotype_metabolite_association/scripts/09_power_analysis.py -- for a grid
of true rho, synthesize a feature correlated with the REAL phenotype residual at that
rho, compute its partial-correlation p-value with the same t-test as
03_corrected_correlation.py, and re-run two-stage BH-FDR against the REAL background
null p-values from that track's full feature scan (not idealized independent uniform
nulls). Fraction of trials clearing Tier1 (|rho_hat|>0.30 & q<0.05) at each true rho
is the empirical power. Each track (cell, supernatant) is powered separately since
they have different n and different background null distributions.
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
    "supernatant": ["species"],
}
RHO_GRID = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
N_TRIALS = int(os.environ.get("N_TRIALS", "300"))
TARGET_POWER = 0.80
RNG_SEED = 20260811
TIER1_RHO_THRESH, TIER1_Q_THRESH = 0.30, 0.05

print("=" * 80)
print("04: POWER ANALYSIS (what rho would Tier1 = |rho|>0.30 & q<0.05 actually detect?)")
print("=" * 80)

rng = np.random.default_rng(RNG_SEED)


def build_design(design_df, covariate_cols):
    Z_parts = [np.ones((len(design_df), 1))]
    for col in covariate_cols:
        dummies = pd.get_dummies(design_df[col].astype(str), drop_first=True)
        Z_parts.append(dummies.values.astype(float))
    Z = np.hstack(Z_parts)
    return Z, np.linalg.matrix_rank(Z), np.linalg.pinv(Z)


def residualize(M, Z, Z_pinv):
    return M - Z @ (Z_pinv @ M)


def power_curve(design_df, covariate_cols, corrected_results_path, label):
    n = len(design_df)
    Z, rank_Z, Z_pinv = build_design(design_df, covariate_cols)
    df = n - rank_Z - 1
    print(f"\n[{label}] n={n}, design rank={rank_Z}, residual df={df}")

    bg = pd.read_csv(corrected_results_path)
    # ANALYSIS_OK[sample-filter]: drops rows with no p-value (constant-feature rows
    # already excluded upstream in 03_corrected_correlation.py never reach this
    # file); asserted non-trivial size immediately below.
    bg_pvals = bg["pval_corrected"].dropna().values
    n_bg = len(bg_pvals)
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert n_bg > 1000, f"[{label}] expected thousands of background null p-values, got {n_bg}"

    y = design_df["mean_auc_rate"].values.astype(float)
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against the same
    # covariates as 03_corrected_correlation.py.
    y_resid = residualize(rankdata(y).reshape(-1, 1), Z, Z_pinv).ravel()
    y_std = (y_resid - y_resid.mean()) / y_resid.std()

    rows = []
    for rho_true in RHO_GRID:
        n_detected = 0
        for _ in range(N_TRIALS):
            noise = rng.standard_normal(n)
            x = rho_true * y_std + np.sqrt(1 - rho_true ** 2) * noise
            x_c = x - x.mean()
            y_c = y_std - y_std.mean()
            r_hat = (x_c @ y_c) / (np.sqrt((x_c ** 2).sum()) * np.sqrt((y_c ** 2).sum()))
            t_stat = r_hat * np.sqrt(df / max(1 - r_hat ** 2, 1e-12))
            p_hat = 2 * t_dist.sf(abs(t_stat), df=df)

            combined_p = np.concatenate([bg_pvals, [p_hat]])
            _, q_combined, _, _ = multipletests(combined_p, method="fdr_bh", alpha=0.05)
            q_synthetic = q_combined[-1]

            if abs(r_hat) > TIER1_RHO_THRESH and q_synthetic < TIER1_Q_THRESH:
                n_detected += 1

        power = n_detected / N_TRIALS
        rows.append({"track": label, "rho_true": rho_true, "power": power,
                     "n_trials": N_TRIALS, "n_background_nulls": int(n_bg)})
    print(f"  power@rho=0.20 -> {[r['power'] for r in rows if r['rho_true'] == 0.20][0]:.2f}, "
          f"power@rho=0.30 -> {[r['power'] for r in rows if r['rho_true'] == 0.30][0]:.2f}")
    return pd.DataFrame(rows)


def min_detectable_rho(power_df, target=TARGET_POWER):
    sub = power_df.sort_values("rho_true")
    rhos, powers = sub["rho_true"].values, sub["power"].values
    if powers[-1] < target:
        return None
    idx = np.searchsorted(powers, target)
    if idx == 0:
        return float(rhos[0])
    r0, r1 = rhos[idx - 1], rhos[idx]
    p0, p1 = powers[idx - 1], powers[idx]
    return float(r0 + (target - p0) * (r1 - r0) / (p1 - p0)) if p1 > p0 else float(r1)


all_power = []
min_rhos = {}
for track, covariate_cols in TRACK_COVARIATES.items():
    design = pd.read_csv(os.path.join(OUT_DIR, track, "sample_design.csv"))
    power_df = power_curve(
        design, covariate_cols,
        os.path.join(OUT_DIR, track, "corrected_all_correlations.csv.gz"), track,
    )
    all_power.append(power_df)
    min_rhos[track] = min_detectable_rho(power_df)

power_all = pd.concat(all_power, ignore_index=True)
power_all.to_csv(os.path.join(OUT_DIR, "power_analysis_curve.csv"), index=False)
print(f"\n✓ wrote {OUT_DIR}/power_analysis_curve.csv")

print(f"\nMinimum detectable rho at {TARGET_POWER:.0%} power:")
for track, r in min_rhos.items():
    print(f"  {track} (n={len(pd.read_csv(os.path.join(OUT_DIR, track, 'sample_design.csv')))}): {r}")

summary = {
    "rho_grid": RHO_GRID,
    "n_trials_per_point": N_TRIALS,
    "tier1_rho_threshold": TIER1_RHO_THRESH,
    "tier1_q_threshold": TIER1_Q_THRESH,
    "target_power": TARGET_POWER,
    "min_detectable_rho": min_rhos,
    "rng_seed": RNG_SEED,
}
with open(os.path.join(OUT_DIR, "power_analysis_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n" + json.dumps(summary, indent=2))

for track, r in min_rhos.items():
    register_value(f"min_detectable_rho_{track}", r, provenance="outputs/power_analysis_summary.json")
