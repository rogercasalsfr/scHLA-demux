#!/bin/bash
#SBATCH --job-name=scHLA-demux
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
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
#   ./run.sh --pools Pool1,Pool2   # only these pools (e.g. one SLURM job per group of pools)
#   ./run.sh --unlock              # remove the lock of a stopped run by hand (see below)
#   sbatch run.sh                  # SLURM: one job with 16 CPUs and 64 GB (Cell Ranger's minimum)
#   sbatch --cpus-per-task=32 --mem=128G run.sh   # a bigger job
#   sbatch --exclusive --mem=0 run.sh             # a whole node, all its CPUs and memory
#   sbatch --partition=long --time=2-00:00:00 run.sh
#   sbatch run.sh --configfile config/other_config.yaml
#
# sbatch options override the ones above. Other arguments are passed to Snakemake. Logs go to
# logs/scHLA-demux_<date_time>.log on a computer, logs/scHLA-demux_<jobid>.out/.err under SLURM.
# (run_slurm.sh is the old name.)

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
# CPUs on this node first: with sbatch --exclusive that is the whole node, not --cpus-per-task.
CORES="${SLURM_CPUS_ON_NODE:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"

if [ ! -f "${REPO_DIR}/workflow/Snakefile" ]; then
  echo "Submit this script from the scHLA-demux folder (no workflow/Snakefile in ${REPO_DIR})" >&2
  exit 1
fi
if [ ! -x "${SNAKEMAKE}" ]; then
  echo "Snakemake not found at ${SNAKEMAKE}; run ./setup.sh first (or set SNAKEMAKE=/path/to/snakemake)" >&2
  exit 1
fi

dry_run=false
unlock=false
args=()   # passed on to Snakemake
pools=()  # --pools Pool1,Pool2, as Snakemake --config (last: --config takes every value after it)
previous=""
for arg in "$@"; do
  case "${previous}" in
    -c|--cores|-j|--jobs) CORES="${arg}" ;;
    --pools) pools=(--config "pools=${arg}") ;;
  esac
  case "${arg}" in
    -n|--dry-run|--dryrun) dry_run=true ;;
    --unlock) unlock=true ;;
    --cores=*) CORES="${arg#--cores=}" ;;
    --pools=*) pools=(--config "pools=${arg#--pools=}") ;;
  esac
  case "${arg}" in
    --pools|--pools=*) ;;
    *) [ "${previous}" = "--pools" ] || args+=("${arg}") ;;
  esac
  previous="${arg}"
done

cd "${REPO_DIR}"

# A stopped run (scancel, time limit, a closed terminal) leaves Snakemake's lock behind, and the next
# run would stop with "Directory cannot be locked". Each run leaves a note in .snakemake/runs/ while
# it runs; when none of those runs is still going, the lock is a leftover and is removed.
RUNS=.snakemake/runs
still_running() {  # note file: host pid slurm_job_id
  local host pid job
  read -r host pid job < "$1" || return 1
  if [ -n "${job}" ] && command -v squeue > /dev/null; then
    [ -n "$(squeue -h -j "${job}" -t PD,R,S,CG,CF -o %i 2> /dev/null)" ]
  elif [ "${host}" = "$(hostname)" ]; then
    kill -0 "${pid}" 2> /dev/null
  else
    return 0  # another computer without SLURM: no way to check, so take it as running
  fi
}
if ! ${unlock} && [ -n "$(ls -A .snakemake/locks 2> /dev/null)" ]; then
  running=""
  for note in "${RUNS}"/*; do
    [ -f "${note}" ] || continue
    if still_running "${note}"; then running+=" ${note##*/}"; else rm -f "${note}"; fi
  done
  # Snakemake started without run.sh (the program, not e.g. tail .snakemake/log/...)
  for pid in $(pgrep -u "$(id -u)" -f '(^|[/ ])snakemake( |$)' || true); do
    [ "$(readlink "/proc/${pid}/cwd" 2> /dev/null)" = "${REPO_DIR}" ] && running+=" pid_${pid}"
  done
  if [ -z "${running}" ]; then
    rm -f .snakemake/locks/*
    echo "Removed the lock left by a stopped run (no other run from this folder is going on)"
  else
    echo "Other runs from this folder are going on (${running# }): they keep their lock. This run" \
         "only stops if it would write the same files (then wait, or use --pools to run other pools)."
  fi
fi

# Not a SLURM job: keep a log like SLURM does, and say how to keep the run going.
if [ -z "${SLURM_JOB_ID:-}" ] && ! ${dry_run} && ! ${unlock}; then
  mkdir -p logs
  LOG="logs/scHLA-demux_$(date '+%Y%m%d_%H%M%S').log"
  exec > >(tee -a "${LOG}") 2>&1
  echo "Running on this computer (no SLURM). Log: ${REPO_DIR}/${LOG}"
  echo "The run stops if this terminal closes: start it in tmux or screen, or as"
  echo "  nohup ./run.sh > /dev/null 2>&1 &     (then follow the log with: tail -f ${LOG})"
  echo "Stopped runs continue where they left off when started again."
fi
if ! ${dry_run} && ! ${unlock}; then
  mkdir -p "${RUNS}"
  NOTE="${RUNS}/$(hostname)_$$"
  echo "$(hostname) $$ ${SLURM_JOB_ID:-}" > "${NOTE}"
  trap 'rm -f "${NOTE}"' EXIT
  trap 'exit 143' TERM INT HUP
fi
echo "[$(date '+%F %T')] scHLA-demux on $(hostname) with ${CORES} cores"

"${SNAKEMAKE}" \
  -s workflow/Snakefile \
  --cores "${CORES}" \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
  "${args[@]}" "${pools[@]}"

echo "[$(date '+%F %T')] scHLA-demux finished"
