#!/usr/bin/env python3
"""Build the cell (C_*) and supernatant (SUP_*) sample/feature tables for the
copper-AUC metabolite association analysis, as two independent, parallel tracks
(not pooled with a sample_type covariate, unlike
analysis/phenotype_metabolite_association/) -- copper resistance could plausibly be
driven by either intracellular defense metabolites or secreted compounds, and the
user asked for both to be tested in parallel from the start (see
.living/decisions.md, "Copper-AUC metabolite association will test cell and
supernatant fractions as parallel tracks, not sequentially").

Data-quality handling (see data/metadata/rhodotorula_auc_copper/provenance.md and
.living/decisions.md for the ingestion-time decisions this implements):

1. `SAMPLE_NAME = "TFCN_17-332M-1"` is duplicated in the copper file, both rows
   mapping to `MS2_SAMPLE_Cell/Supernatant = C_190/SUP_190` but with two different
   `mean_auc_rate` values and two different `Strain ID` values. Per user directive
   (2026-08-11, "skip reconciliation, drop both"), both rows -- and therefore
   C_190/SUP_190 in both tracks -- are dropped permanently, not just pending.
2. 2 rows have `MS2_SAMPLE_Cell`/`MS2_SAMPLE_Supernatant` = literal "No MS2 Data" --
   no metabolomics counterpart; dropped.
3. Species comes from `ATTRIBUTE_species` in the extended strain-traits metadata,
   NOT the copper file's own `SPECIES` column and NOT phase1's `Species` column --
   both of the latter are known-incomplete/inconsistent for this project (see
   analysis/phenotype_metabolite_association/scripts/01_prepare_data.py docstring
   and .living/learnings.md for the recurring strain-ID-collision pattern).
4. A handful of copper-file sample IDs are not present in `phase1_features_filtered`
   at all (already QC-excluded upstream, e.g. blanks/failed-QC samples in Phase 1) --
   these rows are dropped here with counts logged, not silently lost.

Output per track: outputs/<track>/sample_design.csv, outputs/<track>/features_cleaned.csv.gz
"""
import os
import sys

import pandas as pd

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

FEATURES_PATH = os.path.join(ROOT, "analysis", "phase1_features_filtered.csv.gz")
PHENOTYPE_PATH = os.path.join(ROOT, "analysis", "phase1_phenotype_data.csv.gz")
EXT_META_PATH = os.path.join(
    ROOT, "input_data", "MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz"
)
COPPER_PATH = os.path.join(ROOT, "input_data", "Rhodotorula_AUC_copper.20260811.csv.gz")
STRAIN_KEY = "ATTRIBUTE_ID_1"
KNOWN_CONFLICT_SAMPLE_NAME = "TFCN_17-332M-1"  # see module docstring, point 1

TRACKS = {
    "cell": {"copper_id_col": "MS2_SAMPLE_Cell", "sample_type": "cell_pellet"},
    "supernatant": {"copper_id_col": "MS2_SAMPLE_Supernatant", "sample_type": "supernatant"},
}

print("=" * 80)
print("01: PREPARE DATA (cell + supernatant tracks, copper AUC phenotype)")
print("=" * 80)

features_all = pd.read_csv(FEATURES_PATH)
phenotype_all = pd.read_csv(PHENOTYPE_PATH)
ext_meta = pd.read_csv(EXT_META_PATH, sep="\t")
copper = pd.read_csv(COPPER_PATH)

# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis
# convention (fail loudly on broken assumptions, not silently recover).
assert len(features_all) == len(phenotype_all), (
    f"row count mismatch: features={len(features_all)} phenotype={len(phenotype_all)}"
)
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert (phenotype_all.index == features_all.index).all(), "features/phenotype must share row order"
print(f"phase1 features: {features_all.shape}  phenotype: {phenotype_all.shape}")
print(f"copper AUC rows: {len(copper)}")

ext_meta = ext_meta.copy()
ext_meta["sample_id"] = ext_meta["filename"].astype(str)
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert ext_meta["sample_id"].is_unique, "extended metadata sample_id must be unique"

# --- drop the known TFCN_17-332M-1 conflict from the copper table itself, once,
# before splitting into tracks (point 1) ---
n_copper_before = len(copper)
conflict_mask = copper["SAMPLE_NAME"] == KNOWN_CONFLICT_SAMPLE_NAME
n_conflict = int(conflict_mask.sum())
copper = copper.loc[~conflict_mask].copy()
print(f"\nDropped {n_conflict} copper row(s) for '{KNOWN_CONFLICT_SAMPLE_NAME}' "
      f"(known C_190/SUP_190 conflict, permanently excluded per user directive)")

# --- drop "No MS2 Data" rows (point 2) ---
n_no_ms2 = 0
for col in ("MS2_SAMPLE_Cell", "MS2_SAMPLE_Supernatant"):
    n_no_ms2 += int((copper[col] == "No MS2 Data").sum())
no_ms2_mask = (copper["MS2_SAMPLE_Cell"] == "No MS2 Data") | (copper["MS2_SAMPLE_Supernatant"] == "No MS2 Data")
n_no_ms2_rows = int(no_ms2_mask.sum())
copper = copper.loc[~no_ms2_mask].copy()
print(f"Dropped {n_no_ms2_rows} copper row(s) with 'No MS2 Data' in a sample-ID column")
print(f"copper AUC rows remaining after known-issue exclusions: {len(copper)}")

id_cols = ["sample_id", STRAIN_KEY, "ATTRIBUTE_Source", "ATTRIBUTE_species"]

track_summaries = {}
for track_name, cfg in TRACKS.items():
    print("\n" + "-" * 80)
    print(f"TRACK: {track_name}")
    print("-" * 80)
    out_dir = os.path.join(OUT_DIR, track_name)
    os.makedirs(out_dir, exist_ok=True)

    id_col = cfg["copper_id_col"]
    track_copper = copper[["SAMPLE_NAME", id_col, "mean_auc_rate"]].rename(
        columns={id_col: "sample_id"}
    )
    # ANALYSIS_OK[sample-filter]: a copper strain can legitimately have data for
    # one fraction but not the other (upstream assay coverage gap, not this
    # script's doing); rows for the other fraction's placeholder aren't relevant
    # to this track.
    track_copper = track_copper.dropna(subset=["sample_id"])

    dup = track_copper["sample_id"][track_copper["sample_id"].duplicated()]
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert dup.empty, f"[{track_name}] unexpected duplicate sample_id after known-conflict drop: {dup.tolist()}"

    n_before_phase1_join = len(track_copper)
    phase1_track = phenotype_all.loc[phenotype_all["sample_id"].isin(track_copper["sample_id"])]
    joined = track_copper.merge(
        phase1_track[["sample_id", "Library Plate"]], on="sample_id", how="inner", validate="one_to_one"
    )
    n_missing_in_phase1 = n_before_phase1_join - len(joined)
    print(f"copper rows for this track: {n_before_phase1_join}")
    print(f"  not present in phase1_features_filtered (already QC-excluded upstream): {n_missing_in_phase1}")

    merged_ids = joined[["sample_id"]].merge(
        ext_meta[id_cols], on="sample_id", how="left", validate="one_to_one"
    )
    n_join_missing = merged_ids[STRAIN_KEY].isna().sum()
    if n_join_missing:
        raise AssertionError(
            f"[{track_name}] {n_join_missing} samples failed to match extended metadata by "
            "sample_id -- fix the join before proceeding."
        )
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert merged_ids["ATTRIBUTE_Source"].eq(cfg["sample_type"]).all(), (
        f"[{track_name}] unexpected ATTRIBUTE_Source value(s): "
        f"{sorted(set(merged_ids['ATTRIBUTE_Source']) - {cfg['sample_type']})}"
    )

    design = joined.copy()
    design["strain_id"] = merged_ids[STRAIN_KEY].astype(str)
    design["species"] = merged_ids["ATTRIBUTE_species"]

    n_missing_species = int(design["species"].isna().sum())
    # ANALYSIS_OK[sample-filter]: drops rows missing species/plate/AUC after the
    # merge above; count is captured in n_missing_species and printed immediately
    # below, not silently lost.
    design = design.dropna(subset=["species", "Library Plate", "mean_auc_rate"]).reset_index(drop=True)
    print(f"  dropped for missing species/plate/AUC after merge: {n_missing_species}")

    dup_strain = design["strain_id"][design["strain_id"].duplicated()]
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert dup_strain.empty, (
        f"[{track_name}] strain_id must be unique within a track after cleanup: {dup_strain.tolist()}"
    )

    track_features = features_all.loc[
        phenotype_all["sample_id"].isin(design["sample_id"])
    ].reset_index(drop=True)
    # re-align feature rows to design row order via sample_id (features_all shares
    # index/order with phenotype_all, not with design, so re-derive the mapping
    # explicitly rather than assuming positional alignment)
    sample_to_feat_idx = dict(zip(phenotype_all["sample_id"], range(len(phenotype_all))))
    feat_row_order = [sample_to_feat_idx[sid] for sid in design["sample_id"]]
    track_features = features_all.iloc[feat_row_order].reset_index(drop=True)

    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert len(design) == len(track_features), f"[{track_name}] design/features row mismatch"

    n_species = design["species"].nunique()
    print(f"\n  retained: {len(design)} strains, {n_species} species")
    print(design["species"].value_counts().head(10).to_string())

    design.to_csv(os.path.join(out_dir, "sample_design.csv"), index=False)
    track_features.to_csv(os.path.join(out_dir, "features_cleaned.csv.gz"), index=False, compression="gzip")
    print(f"  wrote {out_dir}/sample_design.csv and features_cleaned.csv.gz")

    track_summaries[track_name] = {
        "n_strains": int(len(design)),
        "n_species": int(n_species),
        "n_dropped_missing_in_phase1": int(n_missing_in_phase1),
        "n_dropped_missing_species_plate_auc": int(n_missing_species),
    }
    register_value(f"n_strains_{track_name}", int(len(design)), provenance="scripts/01_prepare_data.py")
    register_value(f"n_species_{track_name}", int(n_species), provenance="scripts/01_prepare_data.py")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
register_value("n_copper_rows_dropped_known_conflict", n_conflict, provenance=f"SAMPLE_NAME={KNOWN_CONFLICT_SAMPLE_NAME}")
register_value("n_copper_rows_dropped_no_ms2_data", n_no_ms2_rows, provenance="MS2_SAMPLE_Cell/Supernatant == 'No MS2 Data'")
for track_name, s in track_summaries.items():
    print(f"{track_name}: {s}")
