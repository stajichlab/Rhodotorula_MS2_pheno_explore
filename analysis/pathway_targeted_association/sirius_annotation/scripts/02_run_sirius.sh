#!/usr/bin/bash -l
#SBATCH -p epyc -c 4 --mem 32gb --out logs/sirius.%A.log

CPU=4

# Run SIRIUS on the 3 pathway-targeted candidate spectra (torularhodin x2,
# shinorine x1). Login confirmed working 2026-08-11 (`sirius login --show`
# shows an active academic subscription) -- unlike the secreted_products run,
# this one should actually produce CSI:FingerID/CANOPUS results, not just
# formula/fingerprint.
set -euo pipefail
cd "$(dirname "$0")/.."

module load sirius

MGF=outputs/sirius_targets.mgf
OUTDIR=outputs/sirius_project
SIRIUS_HEAP_GB="${SIRIUS_HEAP_GB:-10}"

rm -rf "$OUTDIR"

JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G" sirius --cores $CPU --instance-buffer 1 \
  --input "$MGF" --output "$OUTDIR" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  structure \
  canopus \
  write-summaries

echo "SIRIUS project written to $OUTDIR ; summaries under $OUTDIR/summaries/"
