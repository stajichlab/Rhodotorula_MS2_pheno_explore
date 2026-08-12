#!/usr/bin/env python3
"""Build the SIRIUS target list for the 3 pathway-targeted candidate features
(analysis/pathway_targeted_association/outputs/pathway_target_feature_list.csv),
joined to the raw aligned-feature table for a representative raw mzML file to
extract MS2 from. Adapted from
analysis/secreted_products/sirius_annotation/scripts/00_select_targets.py; the
difference is the fallback-file search here considers BOTH C_*/SUP_* samples
(unlike that script, which is scoped to secreted/supernatant-only features) since
these candidates aren't restricted to either compartment.

Output: outputs/sirius_targets.csv (feature_index, row_id, mz, rt_min, charge,
adduct, adduct_rep_file, rep_file_source, mzml_path, has_ms2), same schema as
the secreted_products sirius_annotation targets file for consistency.
"""
import gzip
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANALYSIS_DIR = HERE.parent
REPO_ROOT = ANALYSIS_DIR.parents[2]
RAW_DATA_ROOT = REPO_ROOT.parent  # Rhodotorula_Metabolites/

TARGETS_IN = REPO_ROOT / "analysis" / "pathway_targeted_association" / "outputs" / "pathway_target_feature_list.csv"
ALIGNED_TABLE = REPO_ROOT / "input_data" / "Rhodotorula_MS2_aligned_features_ms2.csv.gz"
OUT = ANALYSIS_DIR / "outputs" / "sirius_targets.csv"


def resolve_mzml_path(rep_file: str) -> str:
    assert rep_file.endswith(".mzML"), f"unexpected adduct_rep_file format: {rep_file!r}"
    if rep_file.startswith("SUP_"):
        return str(RAW_DATA_ROOT / "ExFab_Supernatant" / rep_file)
    if rep_file.startswith("C_"):
        return str(RAW_DATA_ROOT / "mzML" / rep_file)
    raise AssertionError(f"adduct_rep_file does not match known C_/SUP_ prefixes: {rep_file!r}")


def fallback_rep_file(row: dict) -> str:
    """MZmine left adduct_rep_file blank for all 3 pathway-targeted candidates
    (checked directly, not assumed). Fall back to whichever C_*/SUP_* sample has
    the highest peak area for this feature -- these candidates aren't restricted
    to one compartment, unlike secreted_products' SUP_*-only fallback."""
    best_file, best_area = None, -1.0
    for col, val in row.items():
        if not (col.startswith("C_") or col.startswith("SUP_")) or not col.endswith(".mzML Peak area"):
            continue
        if not val:
            continue
        area = float(val)
        if area > best_area:
            best_area = area
            best_file = col[: -len(" Peak area")]
    assert best_file is not None, f"no C_*/SUP_*.mzML peak area found for row ID {row['row ID']}"
    return best_file


def main() -> None:
    with open(TARGETS_IN, newline="") as fh:
        target_ids = [row["raw_position"] for row in csv.DictReader(fh)]
    assert len(target_ids) == 3, f"expected 3 pathway-targeted candidates, got {len(target_ids)}"
    target_id_set = set(target_ids)
    print(f"[targets] {len(target_ids)} target raw_position values loaded from {TARGETS_IN}", file=sys.stderr)

    found = {}
    with gzip.open(ALIGNED_TABLE, "rt", newline="") as fh:
        reader = csv.DictReader(fh)
        for i, row in enumerate(reader):
            key = str(i)
            if key in target_id_set:
                found[key] = row

    missing = target_id_set - found.keys()
    assert not missing, f"{len(missing)} target raw_position values not found in aligned table: {sorted(missing)}"

    fieldnames = [
        "feature_index", "row_id", "mz", "rt_min", "charge", "adduct",
        "adduct_rep_file", "rep_file_source", "mzml_path", "has_ms2",
    ]
    n_fallback = 0
    with open(OUT, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for key in target_ids:
            row = found[key]
            rep_file = row["adduct_rep_file"].strip()
            if rep_file:
                rep_file_source = "mzmine_adduct_rep_file"
            else:
                rep_file = fallback_rep_file(row)
                rep_file_source = "max_peak_area_fallback"
                n_fallback += 1
            writer.writerow({
                "feature_index": key,
                "row_id": row["row ID"],
                "mz": row["row m/z"],
                "rt_min": row["row retention time"],
                "charge": row["charge"],
                "adduct": row["adduct"],
                "adduct_rep_file": rep_file,
                "rep_file_source": rep_file_source,
                "mzml_path": resolve_mzml_path(rep_file),
                "has_ms2": row["has_ms2"],
            })

    print(f"[targets] wrote {len(target_ids)} rows -> {OUT} ({n_fallback} via peak-area fallback)", file=sys.stderr)


if __name__ == "__main__":
    main()
