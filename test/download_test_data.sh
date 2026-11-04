#!/usr/bin/env bash
#
# Downloads the test dataset into test/data/ (about 8 GB, checked against ENA's MD5s):
#
#   CellBench 10x mixture of 3 lung cancer cell lines (NCI-H1975, NCI-H2228, HCC827),
#   Tian et al. 2019, Nature Methods (GEO GSE111108, one lane: SRA SRR6782109).
#   10x 3' v2 gene expression, about 900 cells, with the published demuxlet call for
#   every cell barcode (github.com/LuyiTian/sc_mixology) as the ground truth.
#
# Run again at any time: finished files are kept, interrupted downloads resume.

set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${TEST_DIR}/data"
ENA="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR678/009/SRR6782109"
TRUTH_URL="https://raw.githubusercontent.com/LuyiTian/sc_mixology/7d00ffb0bd043467f365a9358d7312f294f13c9f/misc/demuxlet_output/sc_10X.best"

# file to write | source | MD5 (from ENA)
FILES=(
  "CellBench_R1.fastq.gz|${ENA}/SRR6782109_1.fastq.gz|ce03171af5a1716d90cd166ac080c9f8"
  "CellBench_R2.fastq.gz|${ENA}/SRR6782109_2.fastq.gz|b7be8bcd37c7161e869ad3cbf40e7c4f"
)

mkdir -p "${DATA}"
for entry in "${FILES[@]}"; do
  IFS="|" read -r name url md5 <<< "${entry}"
  target="${DATA}/${name}"
  if [ -f "${target}" ]; then
    echo "${name}: already downloaded"
    continue
  fi
  echo "${name}: downloading"
  progress=(--silent --show-error)  # a progress bar only in a terminal, not in a log
  [ -t 2 ] && progress=(--progress-bar)
  curl -fL --retry 3 -C - "${progress[@]}" -o "${target}.part" "${url}"
  echo "${name}: checking MD5"
  if [ "$(md5sum < "${target}.part" | cut -d' ' -f1)" != "${md5}" ]; then
    echo "${name}: MD5 does not match ENA's; delete ${target}.part and run this script again" >&2
    exit 1
  fi
  mv "${target}.part" "${target}"
done

if [ ! -s "${DATA}/sc_10X.best" ]; then
  echo "sc_10X.best (ground truth): downloading"
  curl -fsSL -o "${DATA}/sc_10X.best" "${TRUTH_URL}"
fi
echo "Test data ready in ${DATA}"
