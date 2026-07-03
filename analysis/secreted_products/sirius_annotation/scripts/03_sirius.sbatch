#!/usr/bin/env bash
#SBATCH --job-name=sirius_secreted
#SBATCH --account=stajichlab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.log
# Submit from analysis/secreted_products/sirius_annotation/: sbatch scripts/03_sirius.sbatch
#
# Runs SIRIUS as its own dedicated allocation instead of inside a shared
# interactive session. An earlier attempt in a 16G interactive job (job
# 26025365) was repeatedly OOM-killed by the cgroup (confirmed via dmesg)
# because that session's memory was already ~9G consumed by unrelated
# concurrent jobs (conda-build, funannotate, raxml-ng) sharing the same
# cgroup — SIRIUS's own -Xmx cap did not help since the ceiling was shared,
# not per-process. A dedicated batch allocation avoids that contention.
set -euo pipefail
# sbatch copies this script into a spool dir before running it, so $0-based
# paths resolve wrong (confirmed: a prior run tried to write to
# /var/spool/slurmd/outputs). Rely on --chdir given at submission time
# (SLURM sets PWD to it); SLURM_SUBMIT_DIR is the fallback for a plain `sh` run.
cd "${SLURM_SUBMIT_DIR:-$PWD}"

module load sirius

MGF=outputs/sirius_targets.mgf
OUTDIR=outputs/sirius_project
SIRIUS_HEAP_GB="${SIRIUS_HEAP_GB:-24}"

rm -rf "$OUTDIR"

JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G" sirius --cores "${SLURM_CPUS_PER_TASK:-8}" \
  --input "$MGF" --output "$OUTDIR" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  structure \
  canopus \
  write-summaries

echo "SIRIUS project written to $OUTDIR"
