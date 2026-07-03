#!/usr/bin/env bash
# Run SIRIUS on the exported MGF of uniquely-secreted MS2 spectra.
#
# formula/fingerprint/canopus/structure need a login (compute node needs
# outbound internet to unijena web services): `module load sirius && sirius login`
# once interactively before this script's structure/canopus steps will produce
# anything beyond "not logged in" warnings. formula + fingerprint (fragmentation
# tree, molecular formula) work fully offline.
set -euo pipefail
cd "$(dirname "$0")/.."

module load sirius

MGF=outputs/sirius_targets.mgf
OUTDIR=outputs/sirius_project

rm -rf "$OUTDIR"

sirius --input "$MGF" --output "$OUTDIR" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  structure \
  canopus \
  write-summaries

echo "SIRIUS project written to $OUTDIR ; summaries under $OUTDIR/summaries/ (or per-compound dirs for older layouts)"
