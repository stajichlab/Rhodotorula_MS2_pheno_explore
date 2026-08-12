#!/usr/bin/env python3
"""Sanity check, both tracks: is mean_auc_rate itself confounded by species or
Library Plate, before any feature-correlation work? This is the same species-first
check that turned out to matter enormously for the color-phenotype analysis
(pooled correlations there were almost entirely a species confound / Simpson's
paradox -- see .living/decisions.md). If copper AUC is itself strongly
species-structured, the species-partial correction used in 03_corrected_correlation.py
is doing real, necessary work, not a formality.

Method: Kruskal-Wallis across species (non-parametric one-way test, robust to the
heavy species-size imbalance -- R. mucilaginosa is ~200/264 rows); same for
Library Plate. Also reports each species' median AUC and n, since a global
Kruskal-Wallis p-value alone doesn't say which species differ or by how much.
"""
import json
import os

import pandas as pd
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "outputs")
TRACKS = ["cell", "supernatant"]
MIN_GROUP_N = 3  # species/plate groups smaller than this are excluded from the KW test (too few for a meaningful group comparison) but still shown in the per-group table

print("=" * 80)
print("02: CONFOUND CHECK (is mean_auc_rate itself species/plate-structured?)")
print("=" * 80)

summary = {}
for track in TRACKS:
    print("\n" + "-" * 80)
    print(f"TRACK: {track}")
    print("-" * 80)
    design = pd.read_csv(os.path.join(OUT_DIR, track, "sample_design.csv"))
    n = len(design)

    track_summary = {"n": n}
    for covariate in ("species", "Library Plate"):
        groups_all = design.groupby(covariate)["mean_auc_rate"]
        sizes = groups_all.size()
        big_enough = sizes[sizes >= MIN_GROUP_N].index
        groups = [design.loc[design[covariate] == g, "mean_auc_rate"].values for g in big_enough]
        n_excluded_groups = int((sizes < MIN_GROUP_N).sum())

        if len(groups) < 2:
            # Real, not a bug: e.g. phase1_phenotype_data.csv.gz records
            # Library Plate = 1.0 for every single SUP_* sample (295/295, verified
            # against the full table, not just this track's subset) -- an upstream
            # data property (plate metadata was apparently only tracked for the
            # cell-pellet batch), not a join error. A covariate with <2 usable
            # groups is constant for this track and cannot be tested or
            # meaningfully regressed out -- flagged here, handled in
            # 03_corrected_correlation.py by dropping it from that track's design.
            print(f"\n  {covariate}: only {len(groups)} group(s) with n>={MIN_GROUP_N} "
                  f"({sizes.to_dict()}) -- constant/unusable as a covariate for this "
                  "track, skipping Kruskal-Wallis")
            track_summary[covariate] = {
                "kruskal_h": None, "kruskal_p": None,
                "n_groups_tested": int(len(groups)), "n_groups_excluded_small": n_excluded_groups,
                "note": "constant or too few groups for this track -- excluded as a covariate",
            }
            continue
        h_stat, p_val = stats.kruskal(*groups)

        print(f"\n  {covariate}: Kruskal-Wallis H={h_stat:.2f}, p={p_val:.4g} "
              f"({len(groups)} groups with n>={MIN_GROUP_N}, {n_excluded_groups} smaller groups excluded from test)")
        summary_tbl = design.groupby(covariate)["mean_auc_rate"].agg(["size", "median", "std"]).sort_values("median")
        print(summary_tbl.to_string())

        track_summary[covariate] = {
            "kruskal_h": float(h_stat),
            "kruskal_p": float(p_val),
            "n_groups_tested": int(len(groups)),
            "n_groups_excluded_small": n_excluded_groups,
        }

    summary[track] = track_summary

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
def fmt(p):
    if p is None:
        return "n/a (constant covariate for this track)"
    return f"p={p:.4g} ({'CONFOUNDED' if p < 0.05 else 'not significant'})"


for track, s in summary.items():
    print(f"{track}: species {fmt(s['species']['kruskal_p'])}, plate {fmt(s['Library Plate']['kruskal_p'])}")

with open(os.path.join(OUT_DIR, "confound_check_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(f"\n✓ wrote {OUT_DIR}/confound_check_summary.json")
