#!/usr/bin/env bash
# Reproduce the SIRIUS MS2 structural-annotation step end to end.
# Step 1 needs plain python3 (conda miniconda 3.9, pandas/csv only).
# Step 2 needs python3.12 with pyteomics/lxml/psims installed
#   (`/usr/bin/python3.12 -m pip install --user pyteomics lxml psims`) —
#   pyteomics 5.x uses `X | Y` union-type syntax that requires Python >=3.10.
# Step 3 needs `module load sirius` (see scripts/02_run_sirius.sh for the
# login caveat on structure/canopus).
set -euo pipefail
cd "$(dirname "$0")"

python3 scripts/00_select_targets.py
/usr/bin/python3.12 scripts/01_export_mgf.py
bash scripts/02_run_sirius.sh

echo "outputs written to outputs/ (sirius_targets.csv, sirius_targets.mgf, mgf_export_summary.csv, sirius_project/)"
