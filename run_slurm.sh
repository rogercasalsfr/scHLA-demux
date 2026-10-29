#!/bin/bash
#SBATCH --job-name=scHLA-demux
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=42
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/scHLA-demux_%j.out
#SBATCH --error=logs/scHLA-demux_%j.err
#
# Runs scHLA-demux for every pool in config/samples.csv, in one SLURM job.
# Submit it from the scHLA-demux folder:
#
#   sbatch run_slurm.sh                  # run everything
#   sbatch --cpus-per-task=1 --mem=4G --partition=normal1 \
#     run_slurm.sh --dry-run             # only check the inputs and show the plan
#   sbatch run_slurm.sh --configfile config/other_config.yaml
#
# Extra arguments are passed to Snakemake. Logs go to logs/scHLA-demux_<jobid>.*

set -euo pipefail

# Snakemake to use. Override with: SNAKEMAKE=/path/to/snakemake sbatch run_slurm.sh
SNAKEMAKE="${SNAKEMAKE:-/home/rcasals/anaconda3/envs/arcashla_clone/bin/snakemake}"

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CORES="${SLURM_CPUS_PER_TASK:-1}"

if [ ! -f "${REPO_DIR}/workflow/Snakefile" ]; then
  echo "Submit this script from the scHLA-demux folder (no workflow/Snakefile in ${REPO_DIR})" >&2
  exit 1
fi
if [ ! -x "${SNAKEMAKE}" ]; then
  echo "Snakemake not found at ${SNAKEMAKE}; set SNAKEMAKE=/path/to/snakemake" >&2
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
