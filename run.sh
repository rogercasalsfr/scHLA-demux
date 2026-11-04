#!/bin/bash
#SBATCH --job-name=scHLA-demux
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/scHLA-demux_%j.out
#SBATCH --error=logs/scHLA-demux_%j.err
#
# Runs scHLA-demux for every pool in config/samples.csv, on any Linux computer or as a SLURM
# job. From the scHLA-demux folder:
#
#   ./run.sh --dry-run             # only check the inputs and show the plan (seconds)
#   ./run.sh                       # run everything on this computer (all its CPUs)
#   ./run.sh --cores 8             # the same with fewer CPUs, e.g. on a shared computer
#   sbatch run.sh                  # SLURM: run everything as one job on a whole node
#   sbatch run.sh --configfile config/other_config.yaml
#   sbatch --partition=highmem run.sh   # pick a partition (sbatch options override the ones above)
#
# Extra arguments are passed to Snakemake. Logs go to logs/scHLA-demux_<date_time>.log on a
# computer, logs/scHLA-demux_<jobid>.out/.err under SLURM. (run_slurm.sh is the old name.)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# sbatch runs a copy of this script from SLURM's spool folder: then the repo is the folder it was submitted
# from. Run directly (also inside another SLURM job, e.g. srun or salloc), the repo is the script's folder.
if [ -f "${SCRIPT_DIR}/workflow/Snakefile" ]; then
  REPO_DIR="${SCRIPT_DIR}"
else
  REPO_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}"
fi
# Snakemake installed by setup.sh. Override with: SNAKEMAKE=/path/to/snakemake ./run.sh
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

dry_run=false
previous=""
for arg in "$@"; do
  case "${arg}" in
    -n|--dry-run|--dryrun) dry_run=true ;;
    --cores=*) CORES="${arg#--cores=}" ;;
  esac
  case "${previous}" in -c|--cores|-j|--jobs) CORES="${arg}" ;; esac
  previous="${arg}"
done

cd "${REPO_DIR}"
# Not a SLURM job: keep a log like SLURM does, and say how to keep the run going.
if [ -z "${SLURM_JOB_ID:-}" ] && ! ${dry_run}; then
  mkdir -p logs
  LOG="logs/scHLA-demux_$(date '+%Y%m%d_%H%M%S').log"
  exec > >(tee -a "${LOG}") 2>&1
  echo "Running on this computer (no SLURM). Log: ${REPO_DIR}/${LOG}"
  echo "The run stops if this terminal closes: start it in tmux or screen, or as"
  echo "  nohup ./run.sh > /dev/null 2>&1 &     (then follow the log with: tail -f ${LOG})"
  echo "Stopped runs continue where they left off when started again."
fi
echo "[$(date '+%F %T')] scHLA-demux on $(hostname) with ${CORES} cores"

"${SNAKEMAKE}" \
  -s workflow/Snakefile \
  --cores "${CORES}" \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
  "$@"

echo "[$(date '+%F %T')] scHLA-demux finished"
