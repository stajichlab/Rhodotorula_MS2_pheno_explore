#!/usr/bin/env python3
"""Permutation null for the corrected partial correlations (03_corrected_correlation.py),
per track. Simpler than
analysis/phenotype_metabolite_association/scripts/03_permutation_null.py's
strain-block design: each track here already has exactly one row per strain (no
C/SUP pairing within a track -- see 01_prepare_data.py), so a plain permutation of
the phenotype residual across rows is already strain-safe; no block structure needed.

Gives a calibrated, non-parametric p-value for the top nominal candidates in each
track (the parametric t-test in 03 assumes independence, which holds here, but a
permutation check is still standard practice per this repo's robust-analysis
convention and directly comparable to the color-phenotype analysis's own check).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

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
N_PERM = 5000
RNG_SEED = 20260811
TOP_N = 25

print("=" * 80)
print("05: PERMUTATION NULL (per track, plain row permutation)")
print("=" * 80)


def build_design(design_df, covariate_cols):
    Z_parts = [np.ones((len(design_df), 1))]
    for col in covariate_cols:
        dummies = pd.get_dummies(design_df[col].astype(str), drop_first=True)
        Z_parts.append(dummies.values.astype(float))
    return np.hstack(Z_parts)


def residualize(M, Z, Z_pinv):
    return M - Z @ (Z_pinv @ M)


overall_summary = {}
for track, covariate_cols in TRACK_COVARIATES.items():
    print("\n" + "-" * 80)
    print(f"TRACK: {track}")
    print("-" * 80)
    out_dir = os.path.join(OUT_DIR, track)
    rng = np.random.default_rng(RNG_SEED)

    design = pd.read_csv(os.path.join(out_dir, "sample_design.csv"))
    features = pd.read_csv(os.path.join(out_dir, "features_cleaned.csv.gz"))
    corrected = pd.read_csv(os.path.join(out_dir, "corrected_all_correlations.csv.gz"))
    n = len(design)

    Z = build_design(design, covariate_cols)
    Z_pinv = np.linalg.pinv(Z)

    feature_values = features.values.astype(float)
    feature_ranks = np.apply_along_axis(rankdata, 0, feature_values)
    # ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py.
    feature_resid_all = residualize(feature_ranks, Z, Z_pinv)

    top = corrected.sort_values("pval_corrected").head(TOP_N)
    feat_idxs = top["feature_index"].to_numpy()
    print(f"permutation-testing top {len(feat_idxs)} nominal candidates x {N_PERM} perms")

    y = design["mean_auc_rate"].values.astype(float)
    y_rank = rankdata(y)
    # ANALYSIS_OK[batch-correction]: same rank-space residualization as 03_corrected_correlation.py.
    y_resid = residualize(y_rank.reshape(-1, 1), Z, Z_pinv).ravel()

    fr = feature_resid_all[:, feat_idxs]
    fr_c = fr - fr.mean(axis=0)
    fr_norm = np.sqrt((fr_c ** 2).sum(axis=0))
    y_c = y_resid - y_resid.mean()
    obs_rho = (fr_c.T @ y_c) / (fr_norm * np.sqrt((y_c ** 2).sum()))

    null_rho = np.empty((N_PERM, len(feat_idxs)))
    for p in range(N_PERM):
        perm = rng.permutation(n)
        yp = y_resid[perm]
        yp_c = yp - yp.mean()
        null_rho[p, :] = (fr_c.T @ yp_c) / (fr_norm * np.sqrt((yp_c ** 2).sum()))

    perm_p = (np.sum(np.abs(null_rho) >= np.abs(obs_rho), axis=0) + 1) / (N_PERM + 1)

    perm_df = pd.DataFrame({
        "feature_index": feat_idxs,
        "rho_corrected": obs_rho,
        "perm_pval": perm_p,
        "null_rho_2p5": np.percentile(null_rho, 2.5, axis=0),
        "null_rho_97p5": np.percentile(null_rho, 97.5, axis=0),
    }).sort_values("perm_pval")
    perm_df.to_csv(os.path.join(out_dir, "permutation_null_results.csv"), index=False)

    n_significant = int((perm_df["perm_pval"] < 0.05).sum())
    print(f"\n✓ wrote {out_dir}/permutation_null_results.csv")
    print(f"candidates with perm_pval<0.05: {n_significant}/{len(perm_df)}")
    print(perm_df.head(10).to_string(index=False))

    track_summary = {
        "n_permutations": N_PERM,
        "n_candidates_tested": int(len(perm_df)),
        "n_candidates_perm_significant_p05": n_significant,
        "rng_seed": RNG_SEED,
    }
    overall_summary[track] = track_summary
    register_value(f"{track}_n_candidates_perm_significant_p05", n_significant, provenance=f"outputs/{track}/permutation_null_results.csv")

with open(os.path.join(OUT_DIR, "permutation_null_summary.json"), "w") as f:
    json.dump(overall_summary, f, indent=2)
print("\n" + json.dumps(overall_summary, indent=2))
