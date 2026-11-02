#!/usr/bin/env bash
#
# Checks an installation end to end on a small public dataset: 3 cell lines from
# 3 people mixed in one 10x library (about 900 cells), with known HLA types and a
# published donor call for every cell. It is a 3' library, so it checks the whole
# pipeline and the clustering, and that no cluster gets a wrong name (HLA naming
# needs 5' data: some clusters stay unnamed). Run it from any folder after ./setup.sh:
#
#   ./test/run_test.sh
#
# It 1) downloads the test data (about 8 GB, only the first time), 2) runs
# scHLA-demux on it (as a SLURM job when sbatch exists, otherwise on this machine;
# about 2-4 hours), and 3) compares the result with the ground truth. It waits for
# the run to finish, so start it inside tmux/screen or with nohup.
# Results: test/results/ (open test/results/qc_report.html).

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"
CONFIG=test/config/config.yaml
PYTHON="${REPO_DIR}/resources/envs/snakemake/bin/python"
[ -x "${PYTHON}" ] || PYTHON=python3

./test/download_test_data.sh

echo "Checking the inputs..."
./run_slurm.sh --configfile "${CONFIG}" --dry-run --quiet > /dev/null

if command -v sbatch > /dev/null && [ -z "${SCHLA_NO_SLURM:-}" ]; then
  echo "Running scHLA-demux on the test data as a SLURM job (logs in logs/)..."
  sbatch --wait run_slurm.sh --configfile "${CONFIG}"
else
  echo "Running scHLA-demux on the test data on this machine..."
  ./run_slurm.sh --configfile "${CONFIG}"
fi

"${PYTHON}" test/check_test.py
