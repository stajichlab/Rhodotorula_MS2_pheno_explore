#!/usr/bin/env python3
"""Broader companion to 01_match_targets.py: matches the expanded target list
(reference_material/pigment_pathway_targets/pigment_pathway_targets_broad.csv --
37 compounds x 8 adducts, including carotenoid esters, ergosterol pathway, broader
MAA/melanin sets, and flavins) against the full 16,332-row raw table at 30 ppm
(widened from 15 ppm) and matching adduct.

This is a SIRIUS-profiling target list, not a phenotype-correlation one -- unlike
01_match_targets.py, this script does not feed into 02_targeted_correlation.py.
Requested after F-006's narrow 3-candidate search turned out to match none of its
hypothesized identities (see .living/findings/color-phenotype-metabolomics.md F-006
update): the goal here is broader chemical profiling coverage via SIRIUS, not another
phenotype-correlation test.
"""
import os

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

TARGETS_PATH = os.path.join(ROOT, "reference_material", "pigment_pathway_targets", "pigment_pathway_targets_broad.csv")
RAW_TABLE_PATH = os.path.join(ROOT, "input_data", "Rhodotorula_MS2_aligned_features_ms2.csv.gz")
PPM_TOLERANCE = 30

print("=" * 80)
print(f"04: BROAD MATCH ({PPM_TOLERANCE} ppm, expanded compound list) FOR SIRIUS PROFILING")
print("=" * 80)

targets = pd.read_csv(TARGETS_PATH)
print(f"target compound x adduct rows: {len(targets)}")

raw = pd.read_csv(RAW_TABLE_PATH, usecols=["row ID", "row m/z", "row retention time", "adduct", "parent_mass", "has_ms2"])
raw = raw.reset_index(drop=True)
raw["raw_position"] = raw.index
print(f"raw aligned table rows: {len(raw)}")

rows = []
for _, t in targets.iterrows():
    candidates = raw.loc[raw["adduct"] == t["adduct"]].copy()
    if candidates.empty:
        continue
    candidates["ppm_error_mz"] = (candidates["row m/z"] - t["expected_mz"]).abs() / t["expected_mz"] * 1e6
    hits = candidates.loc[candidates["ppm_error_mz"] <= PPM_TOLERANCE]
    for _, h in hits.iterrows():
        parent_ppm = None
        if pd.notna(h["parent_mass"]):
            parent_ppm = abs(h["parent_mass"] - t["monoisotopic_mass"]) / t["monoisotopic_mass"] * 1e6
        rows.append({
            "raw_position": int(h["raw_position"]),
            "row_id": int(h["row ID"]),
            "observed_mz": h["row m/z"],
            "observed_adduct": h["adduct"],
            "observed_rt_min": h["row retention time"],
            "has_ms2": bool(h["has_ms2"]),
            "compound": t["compound"],
            "category": t["category"],
            "confidence": t["confidence"],
            "target_adduct": t["adduct"],
            "expected_mz": t["expected_mz"],
            "ppm_error_mz": round(h["ppm_error_mz"], 2),
            "parent_mass_ppm_error": round(parent_ppm, 2) if parent_ppm is not None else None,
        })

matches = pd.DataFrame(rows).sort_values(["category", "compound", "ppm_error_mz"])
matches.to_csv(os.path.join(OUT_DIR, "pathway_target_matches_broad.csv"), index=False)

n_features_matched = matches["raw_position"].nunique() if len(matches) else 0
n_with_ms2 = matches.loc[matches["has_ms2"], "raw_position"].nunique() if len(matches) else 0
print(f"\n✓ {len(matches)} compound-adduct matches within {PPM_TOLERANCE} ppm")
print(f"✓ {n_features_matched} distinct features matched ({n_with_ms2} with has_ms2=True)")
print("\nBy category:")
print(matches.groupby("category")["raw_position"].nunique().to_string())

dedup = matches.drop_duplicates("raw_position")[
    ["raw_position", "row_id", "observed_mz", "observed_rt_min", "compound", "category", "confidence", "ppm_error_mz"]
]
dedup = dedup.sort_values("ppm_error_mz")
dedup.to_csv(os.path.join(OUT_DIR, "pathway_target_feature_list_broad.csv"), index=False)
print(f"\n✓ wrote {OUT_DIR}/pathway_target_feature_list_broad.csv ({len(dedup)} distinct features)")
print("\nTop 20 tightest matches:")
print(dedup.head(20).to_string(index=False))
