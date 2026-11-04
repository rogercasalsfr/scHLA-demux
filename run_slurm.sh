#!/bin/bash
#SBATCH --job-name=scHLA-demux
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/scHLA-demux_%j.out
#SBATCH --error=logs/scHLA-demux_%j.err
#
# Runs scHLA-demux for every pool in config/samples.csv, in one SLURM job that
# takes a whole node (all its CPUs and memory). Submit it from the scHLA-demux folder:
#
#   ./run_slurm.sh --dry-run             # only check the inputs and show the plan (seconds, no SLURM)
#   sbatch run_slurm.sh                  # run everything
#   sbatch run_slurm.sh --configfile config/other_config.yaml
#   sbatch --partition=highmem run_slurm.sh   # pick a partition (sbatch options override the ones below)
#
# Extra arguments are passed to Snakemake. Logs go to logs/scHLA-demux_<jobid>.*

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# sbatch runs a copy of this script from SLURM's spool folder: then the repo is the folder it was submitted
# from. Run directly (also inside another SLURM job, e.g. srun or salloc), the repo is the script's folder.
if [ -f "${SCRIPT_DIR}/workflow/Snakefile" ]; then
  REPO_DIR="${SCRIPT_DIR}"
else
  REPO_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"
fi
# Snakemake installed by setup.sh. Override with: SNAKEMAKE=/path/to/snakemake sbatch run_slurm.sh
SNAKEMAKE="${SNAKEMAKE:-${REPO_DIR}/resources/envs/snakemake/bin/snakemake}"
CORES="${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}"

if [ ! -f "${REPO_DIR}/workflow/Snakefile" ]; then
  echo "Submit this script from the scHLA-demux folder (no workflow/Snakefile in ${REPO_DIR})" >&2
  exit 1
fi
if [ ! -x "${SNAKEMAKE}" ]; then
  echo "Snakemake not found at ${SNAKEMAKE}; run ./setup.sh first (or set SNAKEMAKE=/path/to/snakemake)" >&2
  exit 1
fi

cd "${REPO_DIR}"
echo "[$(date '+%F %T')] scHLA-demux on $(hostname) with ${CORES} cores"

"${SNAKEMAKE}" \
  -s workflow/Snakefile \
  --cores "${CORES}" \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
  "$@"

echo "[$(date '+%F %T')] scHLA-demux finished"
