#!/usr/bin/env bash
# Run SIRIUS on the exported MGF of uniquely-secreted MS2 spectra.
#
# formula/fingerprint/canopus/structure need a login (compute node needs
# outbound internet to unijena web services): `module load sirius && sirius login`
# once interactively before this script's structure/canopus steps will produce
# anything beyond "not logged in" warnings. formula + fingerprint (fragmentation
# tree, molecular formula) work fully offline.
#
# Memory: SIRIUS's launcher sets -XX:MaxRAMPercentage=65 against the *host's*
# visible RAM, not the cgroup/SLURM allocation. On a shared node with a small
# --mem allocation (e.g. an interactive 16G job on a 503G host) this
# under-caps nothing and the JVM heap grows past the cgroup limit -> silent
# OOM-kill partway through (confirmed via dmesg: oom-kill ... task=java).
# Force an explicit, allocation-aware cap with JAVA_OPTS (last -Xmx wins).
set -euo pipefail
cd "$(dirname "$0")/.."

module load sirius

MGF=outputs/sirius_targets.mgf
OUTDIR=outputs/sirius_project
SIRIUS_HEAP_GB="${SIRIUS_HEAP_GB:-10}"  # override for larger SLURM allocations

rm -rf "$OUTDIR"

JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G" sirius --cores 2 --instance-buffer 1 \
  --input "$MGF" --output "$OUTDIR" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  structure \
  canopus \
  write-summaries

echo "SIRIUS project written to $OUTDIR ; summaries under $OUTDIR/summaries/ (or per-compound dirs for older layouts)"
