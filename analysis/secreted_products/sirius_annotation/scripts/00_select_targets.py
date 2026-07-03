#!/usr/bin/env python3
"""Build the SIRIUS target list: the 130 uniquely-secreted features, joined to
the raw aligned-feature table for the MZmine-assigned representative raw file
(``adduct_rep_file``) that MS2 spectra should be extracted from.

Inputs:
  - analysis/secreted_products/outputs/uniquely_secreted_features.csv (130 rows)
  - input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz (16,332 rows, 658 cols)

Output:
  - outputs/sirius_targets.csv: one row per target feature with row ID, m/z,
    RT, charge, adduct, adduct_rep_file, and the raw mzML path to search.
"""
import gzip
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANALYSIS_DIR = HERE.parent
REPO_ROOT = ANALYSIS_DIR.parents[2]
RAW_DATA_ROOT = REPO_ROOT.parent  # Rhodotorula_Metabolites/

TARGETS_IN = REPO_ROOT / "analysis" / "secreted_products" / "outputs" / "uniquely_secreted_features.csv"
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
    """MZmine leaves adduct_rep_file blank for some rows (isotope/adduct-derived
    rows without their own representative-file assignment). Fall back to the
    SUP_* sample with the highest peak area for this feature, since these are
    all *secreted* features and the analysis is scoped to supernatant MS2."""
    best_file, best_area = None, -1.0
    for col, val in row.items():
        if not col.startswith("SUP_") or not col.endswith(".mzML Peak area"):
            continue
        if not val:
            continue
        area = float(val)
        if area > best_area:
            best_area = area
            best_file = col[: -len(" Peak area")]
    assert best_file is not None, f"no SUP_*.mzML peak area found for row ID {row['row ID']}"
    return best_file


def main() -> None:
    with open(TARGETS_IN, newline="") as fh:
        target_ids = [row["feature_index"] for row in csv.DictReader(fh)]
    assert len(target_ids) == 130, f"expected 130 uniquely-secreted targets, got {len(target_ids)}"
    target_id_set = set(target_ids)
    print(f"[targets] {len(target_ids)} target feature_index values loaded from {TARGETS_IN}", file=sys.stderr)

    # feature_index in the secretion-analysis output is the 0-based row position
    # in the aligned table, NOT the MZmine "row ID" column (checked: row 195 ->
    # 'row ID' 205, row 3043 -> 'row ID' 3663 in uniquely_secreted_features.csv,
    # confirming feature_index is positional).
    found = {}
    with gzip.open(ALIGNED_TABLE, "rt", newline="") as fh:
        reader = csv.DictReader(fh)
        for i, row in enumerate(reader):
            key = str(i)
            if key in target_id_set:
                found[key] = row

    missing = target_id_set - found.keys()
    assert not missing, f"{len(missing)} target feature_index values not found in aligned table: {sorted(missing)[:5]}..."

    fieldnames = [
        "feature_index", "row_id", "mz", "rt_min", "charge", "adduct",
        "adduct_rep_file", "rep_file_source", "mzml_path", "has_ms2",
    ]
    n_written = 0
    n_fallback = 0
    with open(OUT, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
        writer.writeheader()
        for feature_index in target_ids:
            row = found[feature_index]
            assert row["has_ms2"].strip().lower() == "true", (
                f"feature_index {feature_index} (row ID {row['row ID']}) has_ms2={row['has_ms2']!r}, "
                "expected true for a uniquely-secreted target"
            )
            rep_file = row["adduct_rep_file"].strip()
            rep_file_source = "mzmine_adduct_rep_file"
            if not rep_file:
                rep_file = fallback_rep_file(row)
                rep_file_source = "max_sup_peak_area_fallback"
                n_fallback += 1
            writer.writerow({
                "feature_index": feature_index,
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
            n_written += 1
    print(f"[targets] {n_fallback}/{n_written} rows used max-peak-area fallback (blank adduct_rep_file)", file=sys.stderr)

    print(f"[targets] wrote {n_written} rows to {OUT}", file=sys.stderr)
    assert n_written == 130


if __name__ == "__main__":
    main()
