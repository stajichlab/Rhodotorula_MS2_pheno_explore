#!/usr/bin/env python3
"""Merge strain/species/sample-type identity onto the Phase 1 feature + phenotype
tables, so downstream steps can properly account for the species confound and the
C/SUP within-strain pairing that Phase 2 (scripts/phase2_correlation_analysis.py)
did not control for. See ../PHENOTYPE_METABOLITE_ASSOCIATION.md for rationale.

Two data-quality problems were found while building this (both real, not script
bugs) and are handled here explicitly rather than silently:

1. The `Species` column in analysis/phase1_phenotype_data.csv.gz (used throughout
   Phase 0-3) is NaN for ~321/590 samples -- almost every SUP_* row, because it was
   populated only from the cell-pellet metadata. This means Phase 0's species-effect
   F-test and any species-stratified analysis were built on a column that was blank
   for most supernatant samples. `ATTRIBUTE_species` in the extended strain-traits
   metadata is populated for both C_* and SUP_* (only 30/590 NaN) and is used here
   instead.
2. One strain ID, "17-332Y-1", maps to two distinct C/SUP filename pairs
   (C_165/SUP_165 and C_269/SUP_269) with SUP_165 and SUP_269 carrying byte-identical
   phenotype values -- an upstream strain-ID collision, not a genuine biological
   replicate. All 4 rows are dropped (documented below) since the strain-block
   permutation design (03_permutation_null.py) requires an unambiguous strain->sample
   mapping.
"""
import os
import sys

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

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
STRAIN_KEY = "ATTRIBUTE_ID_1"
KNOWN_COLLISION_STRAIN = "17-332Y-1"  # see module docstring, point 2

print("=" * 80)
print("01: PREPARE DATA (attach Species / strain / sample-type identity)")
print("=" * 80)

features = pd.read_csv(FEATURES_PATH)
phenotype = pd.read_csv(PHENOTYPE_PATH)
ext_meta = pd.read_csv(EXT_META_PATH, sep="\t")

# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert len(features) == len(phenotype), (
    f"row count mismatch: features={len(features)} phenotype={len(phenotype)}"
)
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert (phenotype.index == features.index).all(), "features/phenotype must share row order"
print(f"features: {features.shape}  phenotype: {phenotype.shape}")

ext_meta = ext_meta.copy()
ext_meta["sample_id"] = ext_meta["filename"].astype(str)
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert ext_meta["sample_id"].is_unique, "extended metadata sample_id must be unique"

id_cols = ["sample_id", STRAIN_KEY, "ATTRIBUTE_Source", "ATTRIBUTE_species"]
# ANALYSIS_OK[join]: one-to-one sample_id join; validate enforced; missing-key
# check immediately below raises rather than silently dropping rows.
merged_ids = phenotype[["sample_id"]].merge(ext_meta[id_cols], on="sample_id", how="left", validate="one_to_one")
n_join_missing = merged_ids[STRAIN_KEY].isna().sum()
if n_join_missing:
    raise AssertionError(
        f"{n_join_missing} of {len(merged_ids)} samples failed to match extended "
        "metadata by sample_id -- fix the join before proceeding."
    )
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert merged_ids["ATTRIBUTE_Source"].isin(["cell_pellet", "supernatant"]).all(), (
    "unexpected ATTRIBUTE_Source value(s): "
    f"{sorted(set(merged_ids['ATTRIBUTE_Source']) - {'cell_pellet', 'supernatant'})}"
)

design = phenotype[
    ["sample_id", "Library Plate", "Median_ColorLab_L*Mean",
     "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
].copy()
design["strain_id"] = merged_ids[STRAIN_KEY].astype(str)
design["sample_type"] = merged_ids["ATTRIBUTE_Source"].map({"cell_pellet": "C", "supernatant": "SUP"})
design["species"] = merged_ids["ATTRIBUTE_species"]  # NOT phenotype['Species'] -- see docstring point 1

n_species_missing_old_col = phenotype["Species"].isna().sum()
n_species_missing_new_col = design["species"].isna().sum()
print(f"Species column (old, phase1_phenotype_data): {n_species_missing_old_col}/{len(phenotype)} missing")
print(f"species column (new, ATTRIBUTE_species):      {n_species_missing_new_col}/{len(design)} missing")

# --- drop rows with an ambiguous strain->sample mapping (docstring point 2) ---
n_before = len(design)
collision_mask = design["strain_id"] == KNOWN_COLLISION_STRAIN
n_collision = int(collision_mask.sum())
design = design.loc[~collision_mask].copy()
features = features.loc[~collision_mask].copy()
print(f"\nDropped {n_collision} rows for strain '{KNOWN_COLLISION_STRAIN}' (ID collision, ambiguous C/SUP pairing)")

# --- drop rows missing any covariate or outcome needed downstream ---
required = ["species", "Library Plate", "Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
missing_mask = design[required].isna().any(axis=1)
n_missing = int(missing_mask.sum())
design = design.loc[~missing_mask].copy()
features = features.loc[~missing_mask].copy()
print(f"Dropped {n_missing} rows missing species/plate/phenotype (required covariates)")

design = design.reset_index(drop=True)
features = features.reset_index(drop=True)
n_after = len(design)
print(f"\nRetained {n_after}/{n_before} samples ({n_before - n_after} dropped total)")

dup_check = design.groupby(["strain_id", "sample_type"]).size()
# ANALYSIS_OK[threshold]: tier thresholds reproduced from the original scripts/phase2_correlation_analysis.py definition (|rho|>0.30/0.25/0.20, q<0.05/0.05/0.10) for direct comparability.
bad = dup_check[dup_check > 1]
# ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
assert bad.empty, f"strain_id/sample_type should be unique per row after cleanup:\n{bad}"

n_strains = design["strain_id"].nunique()
n_paired = design.groupby("strain_id")["sample_type"].apply(lambda s: set(s) == {"C", "SUP"}).sum()
n_species = design["species"].nunique()
print(f"\nstrains: {n_strains} ({n_paired} with both C and SUP)")
print(f"species: {n_species}")
print(design["sample_type"].value_counts().to_string())
print(design["species"].value_counts().to_string())

design.to_csv(os.path.join(OUT_DIR, "sample_design.csv"), index=False)
features.to_csv(os.path.join(OUT_DIR, "features_cleaned.csv.gz"), index=False, compression="gzip")
print(f"\n✓ wrote {OUT_DIR}/sample_design.csv and features_cleaned.csv.gz")

register_value("n_samples_before_cleanup", int(n_before), provenance="analysis/phase1_phenotype_data.csv.gz")
register_value("n_samples_after_cleanup", int(n_after), provenance="scripts/01_prepare_data.py")
register_value("n_samples_dropped_strain_collision", int(n_collision), provenance="strain 17-332Y-1 ID collision, see module docstring")
register_value("n_samples_dropped_missing_covariate", int(n_missing), provenance="missing species/plate/phenotype after switching to ATTRIBUTE_species")
register_value("n_strains_total", int(n_strains), provenance="input_data/...extended_metadata...:ATTRIBUTE_ID_1")
register_value("n_strains_paired_c_sup", int(n_paired), provenance="strains with both C_* and SUP_* samples present")
register_value("n_species", int(n_species), provenance="ATTRIBUTE_species (replaces sparse phase1 Species column)")
