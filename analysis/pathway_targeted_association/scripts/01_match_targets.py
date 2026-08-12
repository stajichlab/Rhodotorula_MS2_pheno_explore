#!/usr/bin/env python3
"""Match the pigment-pathway target compound list
(reference_material/pigment_pathway_targets/pigment_pathway_targets.csv) against
ALL 16,332 features in the raw aligned table, by m/z (15 ppm) and observed adduct --
deliberately NOT restricted to the 7,341 features that survived Phase 1's generic
prevalence (>=10% of samples)/CV(>=0.1) filters.

That restriction was tried first and is wrong for this use case: matching only the
7,341 survivors found just 1 hit (a low-confidence MAA). Matching the full 16,332
finds 3, including two torularhodin candidates (the diagnostic red Rhodotorula
pigment) at 2.5 ppm and 11.5 ppm, both with MS2 spectra, at RT 5.7/6.4 min (plausibly
late-eluting/hydrophobic, consistent with a real carotenoid on a reversed-phase
column) -- neither of which is in the 7,341-feature survivor set. A genuine
pathway/strain-specific pigment is exactly the kind of feature a blanket
>=10%-of-590-samples prevalence filter would discard (it may be produced by only a
subset of species), so re-including the full raw table is the point of "targeted,"
not an error to correct.

Feature-index bookkeeping (three different indices exist in this project's legacy
pipeline -- documented here so it isn't rediscovered the hard way again):
  1. `row ID` in the raw aligned table (input_data/...aligned_features_ms2.csv.gz) --
     MZmine's own ID, NOT contiguous/positional (max 53031 over 16332 rows).
  2. The CSV column header strings in analysis/phase1_features_filtered.csv.gz (and
     phenotype_metabolite_association's features_cleaned.csv.gz) -- confirmed
     empirically (min=1, max=16331 over 7341 surviving columns out of an original
     16332-row raw table) to be the raw table's 0-based ROW POSITION
     (`raw_table.iloc[int(header)]`), not `row ID`.
  3. `feature_index` as used throughout analysis/phenotype_metabolite_association/
     (02-09) -- the POSITIONAL index (0..7340) among the SURVIVING columns only.
This script works directly in the raw table's 0-based row position (`raw_position`,
same numbering as index 2 above but not restricted to survivors), since most of the
compounds of interest here never made it into the survivor set in the first place.
02_targeted_correlation.py pulls raw peak areas by `raw_position` directly from the
aligned table rather than assuming feature_index/features_cleaned.csv.gz coverage.
"""
import os

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

TARGETS_PATH = os.path.join(ROOT, "reference_material", "pigment_pathway_targets", "pigment_pathway_targets.csv")
RAW_TABLE_PATH = os.path.join(ROOT, "input_data", "Rhodotorula_MS2_aligned_features_ms2.csv.gz")
FEATURES_CLEANED_PATH = os.path.join(
    ROOT, "analysis", "phenotype_metabolite_association", "outputs", "features_cleaned.csv.gz"
)
PPM_TOLERANCE = 15

print("=" * 80)
print("01: MATCH PIGMENT-PATHWAY TARGETS TO ALL 16,332 RAW FEATURES (15 ppm, matching adduct)")
print("=" * 80)

targets = pd.read_csv(TARGETS_PATH)
print(f"target compound x adduct rows: {len(targets)}")

# ANALYSIS_OK[positional-access]: only reading the columns needed for matching;
# `usecols` here is a memory optimization, not a silent row/column drop -- all
# 16,332 raw rows are kept.
raw = pd.read_csv(RAW_TABLE_PATH, usecols=["row ID", "row m/z", "row retention time", "adduct", "parent_mass", "has_ms2"])
raw = raw.reset_index(drop=True)
raw["raw_position"] = raw.index  # 0-based row position, see docstring
print(f"raw aligned table rows: {len(raw)}")

feature_headers = pd.read_csv(FEATURES_CLEANED_PATH, nrows=0).columns.tolist()
survivor_positions = {int(h) for h in feature_headers}
print(f"of these, {len(survivor_positions)} survived Phase 1's generic QC filters "
      f"({len(raw) - len(survivor_positions)} did not and are only reachable via this "
      f"full-table search)")

# --- match by adduct + ppm tolerance on observed row m/z ---
rows = []
for _, t in targets.iterrows():
    candidates = raw.loc[raw["adduct"] == t["adduct"]].copy()
    if candidates.empty:
        continue
    candidates["ppm_error_mz"] = (candidates["row m/z"] - t["expected_mz"]).abs() / t["expected_mz"] * 1e6
    hits = candidates.loc[candidates["ppm_error_mz"] <= PPM_TOLERANCE]
    for _, h in hits.iterrows():
        # cross-check against parent_mass (adduct-independent deconvoluted neutral
        # mass), when MZmine reports one, as an independent corroborating signal
        parent_ppm = None
        if pd.notna(h["parent_mass"]):
            parent_ppm = abs(h["parent_mass"] - t["monoisotopic_mass"]) / t["monoisotopic_mass"] * 1e6
        rows.append({
            "raw_position": int(h["raw_position"]),
            "row_id": int(h["row ID"]),
            "survived_phase1_qc": int(h["raw_position"]) in survivor_positions,
            "observed_mz": h["row m/z"],
            "observed_adduct": h["adduct"],
            "observed_rt_min": h["row retention time"],
            "has_ms2": bool(h["has_ms2"]),
            "compound": t["compound"],
            "category": t["category"],
            "pathway_step": t["pathway_step"],
            "confidence": t["confidence"],
            "target_adduct": t["adduct"],
            "expected_mz": t["expected_mz"],
            "ppm_error_mz": round(h["ppm_error_mz"], 2),
            "parent_mass_ppm_error": round(parent_ppm, 2) if parent_ppm is not None else None,
            "notes": t["notes"],
        })

matches = pd.DataFrame(rows).sort_values(["category", "compound", "ppm_error_mz"])
matches.to_csv(os.path.join(OUT_DIR, "pathway_target_matches.csv"), index=False)

n_features_matched = matches["raw_position"].nunique()
print(f"\n✓ {len(matches)} compound-adduct matches within {PPM_TOLERANCE} ppm")
print(f"✓ {n_features_matched} distinct features matched")
print("\nBy category:")
print(matches.groupby("category")["raw_position"].nunique().to_string())
print("\nBy compound:")
print(matches.groupby("compound")["raw_position"].nunique().to_string())
print(f"\nfeatures with MS2 spectra available: {matches.loc[matches['has_ms2'], 'raw_position'].nunique()}/{n_features_matched}")
n_not_in_survivors = (~matches.drop_duplicates('raw_position')['survived_phase1_qc']).sum()
print(f"features NOT in the Phase-1-filtered 7,341 (only reachable via this full-table search): "
      f"{n_not_in_survivors}/{n_features_matched}")

if n_features_matched == 0:
    print("\nNo matches -- nothing to correlation-test. Check target list / ppm tolerance / RT plausibility.")
else:
    print("\nAll matched features (dedup, for the correlation step):")
    dedup = matches.drop_duplicates("raw_position")[
        ["raw_position", "row_id", "survived_phase1_qc", "observed_mz", "observed_rt_min",
         "compound", "category", "confidence", "ppm_error_mz"]
    ]
    dedup.to_csv(os.path.join(OUT_DIR, "pathway_target_feature_list.csv"), index=False)
    print(dedup.to_string(index=False))
