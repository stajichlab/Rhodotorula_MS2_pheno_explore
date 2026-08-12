#!/usr/bin/env bash
#SBATCH --job-name=sirius_broad_profiling
#SBATCH --account=stajichlab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.log
# Submit from analysis/pathway_targeted_association/sirius_annotation_broad/:
#   sbatch --chdir="$PWD" scripts/03_export_and_run.sbatch
#
# Combines the MGF export (01_export_mgf.py, needs python3.12 for pyteomics'
# `X | Y` union-type syntax) and the SIRIUS run into ONE sbatch job so both steps
# persist on the cluster independent of any interactive session -- the earlier
# attempt to run the export interactively in the terminal was killed when the
# session needed to close, before it had written anything. Walltime budgeted
# generously (10h) since 263 spectra (vs. the earlier 3- and 117-spectra runs)
# means both the mzML-scanning export step and the CSI:FingerID/CANOPUS network
# round-trips take proportionally longer.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

echo "=== [1/2] MGF export ==="
/usr/bin/python3.12 scripts/01_export_mgf.py

echo ""
echo "=== [2/2] SIRIUS run ==="
module load sirius

MGF=outputs/sirius_targets.mgf
OUTDIR=outputs/sirius_project
SIRIUS_HEAP_GB="${SIRIUS_HEAP_GB:-24}"

rm -rf "$OUTDIR"

# Subcommand chain fixed 2026-08-11 for SIRIUS 6.x: `structures` now depends on
# CANOPUS output and must come AFTER `canopus`; the 5.x order (...fingerprint,
# structure, canopus...) throws a picocli UnmatchedArgumentException in 6.x.
# -XX:AOTMode=off avoids the shared install-dir AOT cache SIGILL-crashing on
# some of this cluster's CPUs (see .living/learnings.md for both issues).
JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G -XX:AOTMode=off" sirius --cores "${SLURM_CPUS_PER_TASK:-8}" \
  --input "$MGF" --output "$OUTDIR" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  canopus \
  structures \
  write-summaries

echo "SIRIUS project written to $OUTDIR"
