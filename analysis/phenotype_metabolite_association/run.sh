#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

python3 scripts/01_prepare_data.py
python3 scripts/02_corrected_correlation.py
python3 scripts/03_permutation_null.py
python3 scripts/04_replication_holdout.py
python3 scripts/05_within_species_mucilaginosa.py
N_PERM=300 python3 scripts/06_multivariate_module_test.py
python3 scripts/07_mucilaginosa_holdout.py
N_PERM=200 python3 scripts/08_random_forest_module_test.py
N_TRIALS=300 python3 scripts/09_power_analysis.py
