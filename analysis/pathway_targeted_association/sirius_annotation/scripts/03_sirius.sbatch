#!/usr/bin/env bash
#SBATCH --job-name=sirius_pathway_targeted
#SBATCH --account=stajichlab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.log
# Submit from analysis/pathway_targeted_association/sirius_annotation/:
#   sbatch --chdir="$PWD" scripts/03_sirius.sbatch
# Only 3 spectra -- should finish in minutes, but run as a dedicated batch job
# anyway (not interactively) per this project's established practice after the
# OOM/cgroup-sharing issue documented for the secreted_products SIRIUS run.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

module load sirius

MGF=outputs/sirius_targets.mgf
OUTDIR=outputs/sirius_project
SIRIUS_HEAP_GB="${SIRIUS_HEAP_GB:-10}"

rm -rf "$OUTDIR"

JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G" sirius --cores "${SLURM_CPUS_PER_TASK:-4}" \
  --input "$MGF" --output "$OUTDIR" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  structure \
  canopus \
  write-summaries

echo "SIRIUS project written to $OUTDIR"
