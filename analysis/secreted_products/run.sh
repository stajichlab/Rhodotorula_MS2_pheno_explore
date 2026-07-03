#!/usr/bin/env bash
# Reproduce the secreted-products analysis end to end.
# Requires: python3 (conda miniconda 3.9 works) with pandas, numpy, scipy,
# statsmodels, matplotlib. scilintr (optional lint) needs python >=3.11.
set -euo pipefail
cd "$(dirname "$0")"

python3 scripts/01_secretion_analysis.py

echo "outputs written to outputs/ ; figures in outputs/figures/"
