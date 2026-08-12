#!/usr/bin/env bash
#SBATCH --job-name=sirius_secreted
#SBATCH --account=stajichlab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.log
# Submit from analysis/secreted_products/sirius_annotation/: sbatch scripts/03_sirius.sbatch
#
# --time bumped 8h -> 20h 2026-08-11: the two prior attempts (a) failed at
# every CSI:FingerID/CANOPUS step with "Service URL is null" (missing `sirius
# login`, now fixed and confirmed via `sirius login --show`) and (b) the
# second attempt was separately cancelled by Slurm at the 8h wall clock ~3h
# into the CANOPUS step, having only reached 86/131 instance directories.
# Login is fixed, but structure/CANOPUS calls are slow network round-trips
# per compound, so walltime is bumped generously rather than guessing again.
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

# Subcommand chain fixed 2026-08-11 for SIRIUS 6.x: `structures` now depends on
# CANOPUS output and must come AFTER `canopus`; the 5.x order (...fingerprint,
# structure, canopus...) throws a picocli UnmatchedArgumentException in 6.x.
JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G -XX:AOTMode=off" sirius --cores "${SLURM_CPUS_PER_TASK:-8}" \
  --input "$MGF" --output "$OUTDIR" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  canopus \
  structures \
  write-summaries

echo "SIRIUS project written to $OUTDIR"
