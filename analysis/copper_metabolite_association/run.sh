#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

python3 scripts/01_prepare_data.py
python3 scripts/02_confound_check.py
python3 scripts/03_corrected_correlation.py
N_TRIALS=300 python3 scripts/04_power_analysis.py
python3 scripts/05_permutation_null.py
python3 scripts/06_holdout_replication.py
N_PERM=300 python3 scripts/07_multivariate_module_test.py
N_PERM=200 python3 scripts/08_random_forest_module_test.py
# Requires the sibling Rhodotorula_Rodeo repo's
# metal_resistance/results/figures/met_vs_copper_auc_data.csv to already exist (see
# that repo's metal_resistance/scripts/04c_plot_met_correlation.py) -- this run.sh
# does not regenerate that file, only reads it.
python3 scripts/09_cross_reference_genome_signal.py
