#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# Requires analysis/phenotype_metabolite_association/outputs/{sample_design.csv,
# features_cleaned.csv.gz} to already exist (run that analysis's run.sh first).
python3 ../../reference_material/pigment_pathway_targets/build_target_list.py
python3 scripts/01_match_targets.py
python3 scripts/02_targeted_correlation.py
python3 scripts/03_holdout_validation.py
