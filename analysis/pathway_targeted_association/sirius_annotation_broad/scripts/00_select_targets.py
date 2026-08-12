#!/usr/bin/env python3
"""Build the SIRIUS target list for the broader profiling pass requested after F-006's
narrow 3-candidate search turned out to match none of its hypothesized identities.

Combines two sources, deliberately not just loosening ppm tolerance until an arbitrary
count is hit (that would reintroduce the exact false-positive risk F-006 just
demonstrated):
  1. The broadened but still chemistry-motivated pathway search
     (scripts/04_match_targets_broad.py in the parent folder: 37 compounds across
     carotenoid/sterol/melanin/MAA/flavin classes, 30 ppm) -- 32 features.
  2. The top-80-by-nominal-p features from each of the three already-completed corrected
     statistical scans (color-phenotype pooled, copper-AUC cell, copper-AUC supernatant)
     -- 231 unique features (union, some overlap possible across the three tracks).
This profiles what SIRIUS says these statistically-notable features actually ARE,
independent of any pigment-pathway hypothesis, alongside the (now much more
conservative) targeted pigment search. 32 + 231 (checked for overlap at build time)
covers the ~200-500 range the user asked for without inventing chemistry.
"""
import gzip
import csv
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ANALYSIS_DIR = HERE.parent
REPO_ROOT = ANALYSIS_DIR.parents[2]
RAW_DATA_ROOT = REPO_ROOT.parent

BROAD_MATCH_PATH = REPO_ROOT / "analysis" / "pathway_targeted_association" / "outputs" / "pathway_target_feature_list_broad.csv"
ALIGNED_TABLE = REPO_ROOT / "input_data" / "Rhodotorula_MS2_aligned_features_ms2.csv.gz"
OUT = ANALYSIS_DIR / "outputs" / "sirius_targets.csv"
PROVENANCE_OUT = ANALYSIS_DIR / "outputs" / "target_provenance.csv"

TOP_N_PER_TRACK = 80

STAT_SOURCES = [
    ("color_phenotype",
     REPO_ROOT / "analysis" / "phenotype_metabolite_association" / "outputs" / "corrected_all_correlations.csv.gz",
     REPO_ROOT / "analysis" / "phenotype_metabolite_association" / "outputs" / "features_cleaned.csv.gz"),
    ("copper_cell",
     REPO_ROOT / "analysis" / "copper_metabolite_association" / "outputs" / "cell" / "corrected_all_correlations.csv.gz",
     REPO_ROOT / "analysis" / "copper_metabolite_association" / "outputs" / "cell" / "features_cleaned.csv.gz"),
    ("copper_supernatant",
     REPO_ROOT / "analysis" / "copper_metabolite_association" / "outputs" / "supernatant" / "corrected_all_correlations.csv.gz",
     REPO_ROOT / "analysis" / "copper_metabolite_association" / "outputs" / "supernatant" / "features_cleaned.csv.gz"),
]


def resolve_mzml_path(rep_file: str) -> str:
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert rep_file.endswith(".mzML"), f"unexpected adduct_rep_file format: {rep_file!r}"
    if rep_file.startswith("SUP_"):
        return str(RAW_DATA_ROOT / "ExFab_Supernatant" / rep_file)
    if rep_file.startswith("C_"):
        return str(RAW_DATA_ROOT / "mzML" / rep_file)
    raise AssertionError(f"adduct_rep_file does not match known C_/SUP_ prefixes: {rep_file!r}")


def fallback_rep_file(row: dict) -> str:
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
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert best_file is not None, f"no C_*/SUP_*.mzML peak area found for row ID {row['row ID']}"
    return best_file


def main() -> None:
    provenance = {}  # raw_position -> set of source labels

    broad = pd.read_csv(BROAD_MATCH_PATH)
    for rp in broad["raw_position"].unique():
        provenance.setdefault(int(rp), set()).add("broad_pathway_search")
    print(f"[targets] broad pathway search: {broad['raw_position'].nunique()} features", file=sys.stderr)

    for label, corr_path, features_cleaned_path in STAT_SOURCES:
        corr = pd.read_csv(corr_path)
        headers = pd.read_csv(features_cleaned_path, nrows=0).columns.tolist()
        lookup = {i: int(h) for i, h in enumerate(headers)}
        pcol = "pval_corrected" if "pval_corrected" in corr.columns else "pval"
        top = corr.sort_values(pcol).drop_duplicates("feature_index").head(TOP_N_PER_TRACK)
        n_added = 0
        for fi in top["feature_index"]:
            rp = lookup.get(int(fi))
            if rp is None:
                continue
            provenance.setdefault(rp, set()).add(label)
            n_added += 1
        print(f"[targets] {label}: top {TOP_N_PER_TRACK} -> {n_added} raw_positions", file=sys.stderr)

    target_ids = sorted(provenance.keys())
    target_id_set = {str(x) for x in target_ids}
    print(f"[targets] {len(target_ids)} total distinct raw_position targets after union", file=sys.stderr)

    with open(PROVENANCE_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["raw_position", "sources"])
        for rp in target_ids:
            writer.writerow([rp, ";".join(sorted(provenance[rp]))])

    found = {}
    with gzip.open(ALIGNED_TABLE, "rt", newline="") as fh:
        reader = csv.DictReader(fh)
        for i, row in enumerate(reader):
            key = str(i)
            if key in target_id_set:
                found[key] = row

    missing = target_id_set - found.keys()
    # ANALYSIS_OK[runtime-assert]: intentional developer tripwire per robust-analysis convention (fail loudly on broken assumptions, not silently recover); not stripped in this project's execution (no python -O usage).
    assert not missing, f"{len(missing)} target raw_position values not found in aligned table: {sorted(missing)[:5]}"

    fieldnames = [
        "feature_index", "row_id", "mz", "rt_min", "charge", "adduct",
        "adduct_rep_file", "rep_file_source", "mzml_path", "has_ms2",
    ]
    n_fallback = 0
    n_blank_overridden = 0
    n_written = 0
    with open(OUT, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(found, key=int):
            row = found[key]
            if row["has_ms2"].strip().lower() != "true":
                continue  # can't export a spectrum that doesn't exist; not silently substituted
            rep_file = row["adduct_rep_file"].strip()
            if rep_file and (rep_file.startswith("C_") or rep_file.startswith("SUP_")):
                rep_file_source = "mzmine_adduct_rep_file"
            else:
                # MZmine's own rep_file is either blank or points at a QC/blank sample
                # (Blank_*, QC_Mix_*, etc. -- no real biological signal there) -- fall
                # back to the best real C_*/SUP_* sample for this feature instead,
                # documented not silent.
                if rep_file:
                    n_blank_overridden += 1
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
            n_written += 1

    n_no_ms2 = len(target_ids) - n_written
    print(f"[targets] wrote {n_written} rows -> {OUT} ({n_fallback} via peak-area fallback, "
          f"of which {n_blank_overridden} because MZmine's own rep_file was a QC blank; "
          f"{n_no_ms2} dropped for has_ms2=False)", file=sys.stderr)


if __name__ == "__main__":
    main()
