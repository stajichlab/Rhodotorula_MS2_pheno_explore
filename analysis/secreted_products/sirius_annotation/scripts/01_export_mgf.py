#!/usr/bin/env python3
"""Export one representative MS2 spectrum per target feature to MGF for SIRIUS.

For each row in outputs/sirius_targets.csv, scans the assigned mzML file for
MS2 (level-2) spectra whose precursor m/z matches the feature's row m/z within
PPM_TOL and whose scan start time falls within RT_TOL_MIN of the feature's row
retention time. Among matches, keeps the spectrum with the highest total ion
current (most informative fragmentation). If the representative file has no
match, widens the RT window once before giving up (recorded in the summary,
not silently dropped).

Run with the python3.12 interpreter that has pyteomics/lxml/psims installed
(pyteomics 5.x requires the `X | Y` union-type syntax, Python >=3.10).

Output:
  - outputs/sirius_targets.mgf: one entry per successfully matched feature
  - outputs/mgf_export_summary.csv: per-feature match status, for auditing
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

from pyteomics import mzml

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "outputs"
TARGETS_CSV = OUT_DIR / "sirius_targets.csv"
MGF_OUT = OUT_DIR / "sirius_targets.mgf"
SUMMARY_OUT = OUT_DIR / "mgf_export_summary.csv"

PPM_TOL = 15.0
RT_TOL_MIN = 0.15  # +/- minutes, MZmine RT units confirmed == mzML scan start time units
RT_TOL_MIN_WIDE = 0.5
MIN_FRAGMENT_PEAKS = 2  # below this a "spectrum" is noise, not a usable MS2


def ppm_diff(a: float, b: float) -> float:
    return abs(a - b) / b * 1e6


def find_ms2_matches(mzml_path: str, target_mz: float, target_rt: float, rt_tol: float):
    matches = []
    with mzml.MzML(mzml_path) as reader:
        for spec in reader:
            if spec.get("ms level") != 2:
                continue
            scan = spec.get("scanList", {}).get("scan", [{}])[0]
            rt = scan.get("scan start time")
            if rt is None or abs(float(rt) - target_rt) > rt_tol:
                continue
            precursor = spec.get("precursorList", {}).get("precursor", [{}])[0]
            ion = precursor.get("selectedIonList", {}).get("selectedIon", [{}])[0]
            prec_mz = ion.get("selected ion m/z")
            if prec_mz is None or ppm_diff(float(prec_mz), target_mz) > PPM_TOL:
                continue
            mz_arr = spec.get("m/z array")
            int_arr = spec.get("intensity array")
            if mz_arr is None or int_arr is None or len(mz_arr) < MIN_FRAGMENT_PEAKS:
                continue
            matches.append({
                "scan_id": spec.get("id"),
                "rt": float(rt),
                "precursor_mz": float(prec_mz),
                "tic": float(sum(int_arr)),
                "mz_array": list(mz_arr),
                "intensity_array": list(int_arr),
            })
    return matches


def write_mgf_entry(fh, row: dict, spec: dict) -> None:
    fh.write("BEGIN IONS\n")
    fh.write(f"FEATURE_ID={row['row_id']}\n")
    fh.write(f"TITLE=feature_index_{row['feature_index']}_row_{row['row_id']}\n")
    fh.write(f"PEPMASS={spec['precursor_mz']:.6f}\n")
    charge = row["charge"].strip()
    fh.write(f"CHARGE={charge}+\n")
    fh.write(f"RTINSECONDS={spec['rt'] * 60:.2f}\n")
    for mz, inten in zip(spec["mz_array"], spec["intensity_array"]):
        fh.write(f"{mz:.6f} {inten:.2f}\n")
    fh.write("END IONS\n\n")


def main() -> None:
    with open(TARGETS_CSV, newline="") as fh:
        targets = list(csv.DictReader(fh))
    assert len(targets) == 130, f"expected 130 targets, found {len(targets)}"

    summary_rows = []
    n_exported = 0
    with open(MGF_OUT, "w") as mgf_fh:
        for row in targets:
            target_mz = float(row["mz"])
            target_rt = float(row["rt_min"])
            matches = find_ms2_matches(row["mzml_path"], target_mz, target_rt, RT_TOL_MIN)
            widened = False
            if not matches:
                matches = find_ms2_matches(row["mzml_path"], target_mz, target_rt, RT_TOL_MIN_WIDE)
                widened = True
            if not matches:
                summary_rows.append({
                    "feature_index": row["feature_index"], "row_id": row["row_id"],
                    "status": "no_ms2_match", "mzml_path": row["mzml_path"],
                    "n_candidate_spectra": 0, "chosen_scan_id": "", "chosen_tic": "",
                    "rt_window_widened": widened,
                })
                print(f"[export] WARNING no MS2 match for row_id={row['row_id']} in {row['mzml_path']}", file=sys.stderr)
                continue
            best = max(matches, key=lambda m: m["tic"])
            write_mgf_entry(mgf_fh, row, best)
            n_exported += 1
            summary_rows.append({
                "feature_index": row["feature_index"], "row_id": row["row_id"],
                "status": "exported", "mzml_path": row["mzml_path"],
                "n_candidate_spectra": len(matches), "chosen_scan_id": best["scan_id"],
                "chosen_tic": best["tic"], "rt_window_widened": widened,
            })

    with open(SUMMARY_OUT, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"[export] {n_exported}/{len(targets)} features exported to {MGF_OUT}", file=sys.stderr)
    print(f"[export] summary written to {SUMMARY_OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
