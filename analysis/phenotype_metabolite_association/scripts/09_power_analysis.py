#!/usr/bin/env python3
"""Power analysis: what true correlation strength (rho) would this pipeline's primary
detection criterion (partial correlation, Tier1 = |rho_hat|>0.30 AND two-stage BH
q<0.05, exactly as in 02_corrected_correlation.py) actually be able to detect, at the
sample sizes used in this analysis?

This directly answers the open question left after five converging negative checks:
is the null result meaningful, or is the design simply underpowered to see a real but
modest effect?

Method (simulation, not a closed-form formula, so it reflects the true multiple-testing
burden this pipeline imposes): for a grid of true rho values, repeatedly synthesize one
additional "feature" whose rank-residual is correlated with the real, already-computed
phenotype residual at exactly that rho (plus independent noise), compute its partial
correlation p-value with the same t-test used in 02/05, and re-run the two-stage BH-FDR
using the REAL background null p-values from the 7,341-feature scan (not idealized
independent uniform nulls) so the simulated signal has to compete against the actual
observed null distribution -- exactly as any real feature would have. The fraction of
trials where the synthetic feature clears Tier1 is the empirical power at that rho.

Two designs are powered separately, matching the two completed analyses:
  - "pooled": the corrected 550-sample, 17-species model (02_corrected_correlation.py).
  - "mucilaginosa": the 415-sample, single-species model (05_within_species_mucilaginosa.py).

Limitation (documented, not hidden): the synthetic signal is generated directly in
rank-residual space rather than as a raw feature run through the full rank-transform,
which is an adequate approximation for power purposes (rank transform of a synthesized
near-Gaussian signal preserves its correlation structure to good approximation at
n>400) but is not a byte-for-byte replay of the real pipeline. Permutation-null power
(as opposed to the parametric-FDR power computed here) was not simulated -- it would
require re-running ~5,000 permutations per trial per rho value, which is computationally
prohibitive at this grid resolution; the FDR criterion is the one that actually
determined "0/12,269 hits survive," so it is the one powered here.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist
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
RHO_GRID = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
N_TRIALS = int(os.environ.get("N_TRIALS", "300"))
TARGET_POWER = 0.80
RNG_SEED = 20260811
# ANALYSIS_OK[threshold]: Tier1 definition reproduced exactly from
# 02_corrected_correlation.py / scripts/phase2_correlation_analysis.py for direct
# comparability -- this IS the criterion being powered, not an arbitrary choice here.
TIER1_RHO_THRESH = 0.30
TIER1_Q_THRESH = 0.05

print("=" * 80)
print("09: POWER ANALYSIS (what rho would Tier1 = |rho|>0.30 & q<0.05 actually detect?)")
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
    from scipy.stats import rankdata

    n = len(design_df)
    Z, rank_Z, Z_pinv = build_design(design_df, covariate_cols)
    df = n - rank_Z - 1
    print(f"\n[{label}] n={n}, design rank={rank_Z}, residual df={df}")

    bg = pd.read_csv(corrected_results_path)
    bg_col = "pval_corrected" if "pval_corrected" in bg.columns else "pval"

    rows = []
    for phenotype in PHENOTYPE_COLS:
        y = design_df[phenotype].values.astype(float)
        # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against the
        # same covariates as 02_corrected_correlation.py / 05_within_species_mucilaginosa.py;
        # documented in the module docstring.
        y_resid = residualize(rankdata(y).reshape(-1, 1), Z, Z_pinv).ravel()
        y_std = (y_resid - y_resid.mean()) / y_resid.std()

        # ANALYSIS_OK[sample-filter]: selects this phenotype's rows from the
        # background null-correlation table (one row per feature per phenotype);
        # not a data-quality drop, just a per-phenotype slice for the power loop.
        bg_pvals = bg.loc[bg["phenotype"] == phenotype, bg_col].dropna().values
        n_bg = len(bg_pvals)
        # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per
        # robust-analysis convention (fail loudly if the background null table is
        # unexpectedly small/empty); not stripped in this project's execution.
        assert n_bg > 1000, f"expected thousands of background null p-values, got {n_bg}"

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
            rows.append({"design": label, "phenotype": phenotype, "rho_true": rho_true,
                         "power": power, "n_trials": N_TRIALS, "n_background_nulls": int(n_bg)})
        print(f"  {phenotype}: power@rho=0.20 -> "
              f"{[r['power'] for r in rows if r['phenotype'] == phenotype and r['rho_true'] == 0.20][0]:.2f}, "
              f"power@rho=0.30 -> {[r['power'] for r in rows if r['phenotype'] == phenotype and r['rho_true'] == 0.30][0]:.2f}")
    return pd.DataFrame(rows)


def min_detectable_rho(power_df, target=TARGET_POWER):
    """Linear-interpolate the grid to find rho at target power, per phenotype."""
    out = {}
    for phenotype, sub in power_df.groupby("phenotype"):
        sub = sub.sort_values("rho_true")
        rhos, powers = sub["rho_true"].values, sub["power"].values
        if powers[-1] < target:
            out[phenotype] = None  # not reached even at the top of the grid
            continue
        idx = np.searchsorted(powers, target)
        if idx == 0:
            out[phenotype] = float(rhos[0])
        else:
            r0, r1 = rhos[idx - 1], rhos[idx]
            p0, p1 = powers[idx - 1], powers[idx]
            out[phenotype] = float(r0 + (target - p0) * (r1 - r0) / (p1 - p0)) if p1 > p0 else float(r1)
    return out


design_all = pd.read_csv(os.path.join(OUT_DIR, "sample_design.csv"))

print("\n--- Pooled corrected design (550 samples, 17 species) ---")
pooled_power = power_curve(
    design_all, ("species", "Library Plate", "sample_type"),
    os.path.join(OUT_DIR, "corrected_all_correlations.csv.gz"), "pooled",
)

print("\n--- R. mucilaginosa-only design (415 samples) ---")
muc_design = design_all.loc[design_all["species"] == "Rhodotorula mucilaginosa"].reset_index(drop=True)
muc_power = power_curve(
    muc_design, ("Library Plate", "sample_type"),
    os.path.join(OUT_DIR, "mucilaginosa_all_correlations.csv.gz"), "mucilaginosa",
)

power_df = pd.concat([pooled_power, muc_power], ignore_index=True)
power_df.to_csv(os.path.join(OUT_DIR, "power_analysis_curve.csv"), index=False)
print(f"\n✓ wrote {OUT_DIR}/power_analysis_curve.csv")

min_rho_pooled = min_detectable_rho(pooled_power)
min_rho_muc = min_detectable_rho(muc_power)

print(f"\nMinimum detectable rho at {TARGET_POWER:.0%} power:")
print("  pooled (n=550):      ", min_rho_pooled)
print("  mucilaginosa (n=415):", min_rho_muc)

summary = {
    "rho_grid": RHO_GRID,
    "n_trials_per_point": N_TRIALS,
    "tier1_rho_threshold": TIER1_RHO_THRESH,
    "tier1_q_threshold": TIER1_Q_THRESH,
    "target_power": TARGET_POWER,
    "min_detectable_rho_pooled_n550": min_rho_pooled,
    "min_detectable_rho_mucilaginosa_n415": min_rho_muc,
    "rng_seed": RNG_SEED,
}
with open(os.path.join(OUT_DIR, "power_analysis_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n" + json.dumps(summary, indent=2))

for phenotype, rho in min_rho_pooled.items():
    key = phenotype.replace("Median_ColorLab_", "").replace("*Mean", "").lower()
    register_value(f"min_detectable_rho_pooled_{key}", rho, provenance="outputs/power_analysis_summary.json")
for phenotype, rho in min_rho_muc.items():
    key = phenotype.replace("Median_ColorLab_", "").replace("*Mean", "").lower()
    register_value(f"min_detectable_rho_mucilaginosa_{key}", rho, provenance="outputs/power_analysis_summary.json")
