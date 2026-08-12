#!/usr/bin/env bash
#SBATCH --job-name=sirius_broad_profiling
#SBATCH --account=stajichlab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.log
# Submit from analysis/pathway_targeted_association/sirius_annotation/:
#   sbatch --chdir="$PWD" scripts/03_sirius.sbatch
# Only 3 spectra -- should finish in minutes, but run as a dedicated batch job
# anyway (not interactively) per this project's established practice after the
# OOM/cgroup-sharing issue documented for the secreted_products SIRIUS run.
#
# Subcommand chain fixed 2026-08-11 for SIRIUS 6.x (was formula/fingerprint/
# structure/canopus/write-summaries, the 5.x order): `structures` now depends
# on CANOPUS output and must come AFTER `canopus`, not before it, and the
# 5.x-era singular alias `structure` produces a picocli UnmatchedArgumentException
# in 6.x (confirmed via `sirius formula|fingerprint|canopus --help`'s own
# Commands: listings -- the valid chain is formula -> fingerprint -> canopus ->
# structures -> write-summaries).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

module load sirius

MGF=outputs/sirius_targets.mgf
OUTDIR=outputs/sirius_project
SIRIUS_HEAP_GB="${SIRIUS_HEAP_GB:-24}"

rm -rf "$OUTDIR"

JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G -XX:AOTMode=off" sirius --cores "${SLURM_CPUS_PER_TASK:-4}" \
  --input "$MGF" --output "$OUTDIR" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  canopus \
  structures \
  write-summaries

echo "SIRIUS project written to $OUTDIR"
