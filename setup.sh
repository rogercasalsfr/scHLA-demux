#!/usr/bin/env bash
#
# Installs everything scHLA-demux needs into resources/. Run it once after
# downloading the repo, from any folder:
#
#   ./setup.sh             # asks for the Cell Ranger download link, installs the rest
#   ./setup.sh --cellranger-url "<link>"   # the same, without the question
#   ./setup.sh --check     # only report what is installed
#
# Run it again at any time: finished parts are skipped. Everything is public
# except Cell Ranger, which 10x only gives out after you accept their license
# (setup.sh explains how). Works on a cluster and on a single Linux computer:
# without conda, Miniforge is installed into resources/ first, and without
# apptainer (or singularity), apptainer is installed there too, with no admin
# rights. config/config.yaml already points to resources/, so nothing needs
# editing afterwards.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RES="${REPO_DIR}/resources"

GEX_REF="refdata-gex-GRCh38-2024-A"
GEX_URL="https://cf.10xgenomics.com/supp/cell-exp/${GEX_REF}.tar.gz"
VDJ_REF="refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
VDJ_URL="https://cf.10xgenomics.com/supp/cell-vdj/${VDJ_REF}.tar.gz"
# souporcell images, tried in order: the Cumulus (Broad) Docker build of souporcell
# 2.5, then the author's image on the Sylabs library.
SOUPORCELL_IMAGES=(
  "docker://cumulusprod/souporcell:2.5"
  "library://wheaton5/souporcell/souporcell:release"
)
SYLABS_LIBRARY="https://library.sylabs.io"
# IPD-IMGT/HLA release used by arcasHLA: 3.66.0 (tag v3.66.0-alpha), the release the
# pipeline is tested with. HLA allele names and sequences change between releases.
IMGT_RELEASE="3.66.0"
IMGT_COMMIT="d972ae287948dc2e5468a40ff3b3d27227541189"
SUBSET_BAM_URL="https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux"
CELLRANGER_PAGE="https://www.10xgenomics.com/support/software/cell-ranger/downloads"

CELLRANGER_URL="${CELLRANGER_URL:-}"
CHECK_ONLY=false

usage() {
  sed -n '3,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [ $# -gt 0 ]; do
  case "$1" in
    --cellranger-url) CELLRANGER_URL="${2:-}"; shift 2 ;;
    --check) CHECK_ONLY=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

newest() {  # newest existing path among the arguments, skipping unfinished downloads
  local path
  for path in "$@"; do
    case "$path" in *.part) continue ;; esac
    if [ -e "$path" ]; then echo "$path"; fi
  done | sort -V | tail -n 1
}

download() {  # url destination
  log "downloading $(basename "$2")"
  curl -fL --retry 3 -C - -o "$2.part" "$1"
  mv "$2.part" "$2"
}

unpack() {  # archive folder_name: unpacks into resources/<folder_name>, all or nothing
  local archive=$1 name=$2 tmp="${RES}/.unpack_$2"
  log "unpacking $(basename "${archive}")"
  rm -rf "${tmp}"
  mkdir -p "${tmp}"
  if ! tar -xf "${archive}" -C "${tmp}"; then
    rm -rf "${tmp}"
    return 1
  fi
  if [ ! -d "${tmp}/${name}" ]; then
    echo "$(basename "${archive}") does not contain a ${name}/ folder" >&2
    rm -rf "${tmp}"
    return 1
  fi
  mv "${tmp}/${name}" "${RES}/${name}"
  rm -rf "${tmp}"
}

# --------------------------------------------------------------------------
# Cell Ranger
# --------------------------------------------------------------------------

cellranger_help() {
  cat >&2 <<EOF
Cell Ranger needs a download link from 10x Genomics (free, but you have to
accept their license first):
  1. Open ${CELLRANGER_PAGE}
  2. Fill in the form. The page then shows a curl command.
  3. Run ./setup.sh again and paste that command (or just the link in it) when
     it asks. Or give the link on the command line, inside double quotes (the
     link contains '&'): ./setup.sh --cellranger-url "<the link>"
The link expires after about a day and is personal: never put it in a file of
the repo. If you already downloaded cellranger-<version>.tar.gz or .tar.xz,
move it into resources/ and run ./setup.sh again.
EOF
}

ask_cellranger_url() {  # in a terminal: ask for the link, accepting the whole curl command
  local answer
  [ -t 0 ] || return 0
  cat <<EOF

Cell Ranger is free, but 10x Genomics gives the download link only after you
accept their license:
  1. Open ${CELLRANGER_PAGE}
  2. Fill in the form. The page then shows a curl command.
Paste that command (or just the link in it) here and press Enter. Press Enter
alone to install everything else now and Cell Ranger later.
EOF
  read -r -p "> " answer || true
  CELLRANGER_URL=$(grep -o "https://cf\.10xgenomics\.com/releases/[^\"' ]*" <<< "${answer}" | head -n 1 || true)
  if [ -n "${answer}" ] && [ -z "${CELLRANGER_URL}" ]; then
    echo "That does not contain a Cell Ranger link (https://cf.10xgenomics.com/releases/...)." >&2
  fi
}

check_cellranger() {
  [ -x "${RES}/cellranger/cellranger" ] && "${RES}/cellranger/cellranger" --version >/dev/null 2>&1
}

install_cellranger() {
  local exe archive name
  exe=$(newest "${RES}"/cellranger-*/cellranger)
  if [ -z "${exe}" ]; then
    archive=$(newest "${RES}"/cellranger-*.tar.*)
    if [ -z "${archive}" ] && [ -z "${CELLRANGER_URL}" ]; then
      ask_cellranger_url
    fi
    if [ -z "${archive}" ] && [ -n "${CELLRANGER_URL}" ]; then
      archive="${RES}/$(basename "${CELLRANGER_URL%%\?*}")"
      case "$(basename "${archive}")" in
        cellranger-*.tar.*) ;;
        *) echo "This is not a Cell Ranger download link: ${CELLRANGER_URL%%\?*}" >&2; return 1 ;;
      esac
      download "${CELLRANGER_URL}" "${archive}"
    fi
    if [ -z "${archive}" ]; then
      cellranger_help
      return 1
    fi
    name=$(basename "${archive}")
    name=${name%%.tar.*}
    if ! unpack "${archive}" "${name}"; then
      echo "$(basename "${archive}") is damaged: delete it and run setup again with a new link" >&2
      return 1
    fi
    exe="${RES}/${name}/cellranger"
  fi
  ln -sfn "$(basename "$(dirname "${exe}")")" "${RES}/cellranger"
}

# --------------------------------------------------------------------------
# 10x references
# --------------------------------------------------------------------------

gex_ok() {
  [ -f "$1/reference.json" ] && [ -s "$1/fasta/genome.fa" ] && [ -s "$1/star/SA" ] \
    && compgen -G "$1/genes/genes.gtf*" >/dev/null
}

vdj_ok() {
  [ -f "$1/reference.json" ] && [ -s "$1/fasta/regions.fa" ]
}

install_reference() {  # folder_name url check_function link_name
  local name=$1 url=$2 check=$3 link=$4 archive="${RES}/$1.tar.gz"
  if ! "${check}" "${RES}/${name}"; then
    if [ -e "${RES}/${name}" ]; then
      log "${name} is incomplete, unpacking it again"
      rm -rf "${RES:?}/${name}"
    fi
    [ -f "${archive}" ] || download "${url}" "${archive}"
    if ! unpack "${archive}" "${name}"; then
      log "$(basename "${archive}") is damaged, downloading it again"
      rm -f "${archive}"
      download "${url}" "${archive}"
      unpack "${archive}" "${name}"
    fi
  fi
  ln -sfn "${name}" "${RES}/${link}"
}

check_gex() { gex_ok "${RES}/refdata-gex"; }
install_gex() { install_reference "${GEX_REF}" "${GEX_URL}" gex_ok refdata-gex; }
check_vdj() { vdj_ok "${RES}/refdata-vdj"; }
install_vdj() { install_reference "${VDJ_REF}" "${VDJ_URL}" vdj_ok refdata-vdj; }

# --------------------------------------------------------------------------
# souporcell and subset-bam
# --------------------------------------------------------------------------

# apptainer (or singularity) runs the souporcell image. Clusters usually have it; on a computer
# without it, conda-forge's apptainer goes into resources/envs/apptainer (no admin rights needed).
find_apptainer() {  # same order as the workflow: an installed copy first, then the system's
  if [ -x "${RES}/envs/apptainer/bin/apptainer" ]; then
    APPTAINER="${RES}/envs/apptainer/bin/apptainer"
  else
    APPTAINER=$(command -v apptainer || command -v singularity || true)
  fi
}
find_apptainer
SOUPORCELL_SIF="${RES}/souporcell_release.sif"

check_apptainer() {  # looks again: the install runs in a subshell
  find_apptainer
  [ -n "${APPTAINER}" ] && "${APPTAINER}" --version >/dev/null 2>&1
}

install_apptainer() {
  log "apptainer (or singularity) not found: installing apptainer into resources/envs/apptainer"
  install_env apptainer
}

check_souporcell() {
  [ -n "${APPTAINER}" ] && [ -s "${SOUPORCELL_SIF}" ] \
    && "${APPTAINER}" exec "${SOUPORCELL_SIF}" souporcell_pipeline.py -h >/dev/null 2>&1
}

install_souporcell() {
  if [ -z "${APPTAINER}" ]; then
    echo "apptainer is not available (see the apptainer step above)" >&2
    return 1
  fi
  local tmp="${RES}/.souporcell_download.sif" image
  for image in "${SOUPORCELL_IMAGES[@]}"; do
    log "downloading the souporcell image (about 2 GB) from ${image}"
    rm -f "${tmp}"
    if "${APPTAINER}" pull --library "${SYLABS_LIBRARY}" "${tmp}" "${image}" \
      && "${APPTAINER}" exec "${tmp}" souporcell_pipeline.py -h >/dev/null 2>&1; then
      mv "${tmp}" "${SOUPORCELL_SIF}"
      return 0
    fi
    log "could not use ${image}"
  done
  rm -f "${tmp}"
  if [ "$(cat /proc/sys/kernel/apparmor_restrict_unprivileged_userns 2>/dev/null)" = 1 ]; then
    echo "This system blocks the user namespaces that apptainer needs without admin rights (Ubuntu 23.10" >&2
    echo "and later). Ask an administrator to install apptainer system-wide (its package allows it), then" >&2
    echo "run ./setup.sh again." >&2
  fi
  return 1
}

check_subset_bam() {
  [ -x "${RES}/bin/subset-bam" ] && "${RES}/bin/subset-bam" --version >/dev/null 2>&1
}

install_subset_bam() {
  mkdir -p "${RES}/bin"
  download "${SUBSET_BAM_URL}" "${RES}/bin/subset-bam"
  chmod +x "${RES}/bin/subset-bam"
}

# --------------------------------------------------------------------------
# conda envs: snakemake, and arcasHLA + samtools
# --------------------------------------------------------------------------

MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
CONDA=$(command -v mamba || command -v conda || true)
if [ -z "${CONDA}" ] && [ -x "${RES}/miniforge/bin/mamba" ]; then
  CONDA="${RES}/miniforge/bin/mamba"
fi

install_miniforge() {  # no conda on this machine: conda-forge's Miniforge into resources/miniforge
  if [ ! -x "${RES}/miniforge/bin/mamba" ]; then
    log "conda not found: installing Miniforge into resources/miniforge (no admin rights needed)"
    download "${MINIFORGE_URL}" "${RES}/miniforge_installer.sh"
    rm -rf "${RES}/miniforge"
    bash "${RES}/miniforge_installer.sh" -b -p "${RES}/miniforge"
    rm -f "${RES}/miniforge_installer.sh"
  fi
  CONDA="${RES}/miniforge/bin/mamba"
}

env_ready() {  # name: installed from the current envs/<name>.yaml
  [ "$(cat "${RES}/envs/$1/.schla_ready" 2>/dev/null)" = "$(md5sum < "${REPO_DIR}/envs/$1.yaml")" ]
}

install_env() {  # name (envs/<name>.yaml -> resources/envs/<name>)
  local prefix="${RES}/envs/$1"
  env_ready "$1" && return 0
  [ -n "${CONDA}" ] || install_miniforge
  rm -rf "${prefix}"
  log "creating the $1 conda env (a few minutes)"
  # no "Confirm changes [Y/n]" question from mamba in the middle of the setup
  CONDA_ALWAYS_YES=true MAMBA_ALWAYS_YES=true \
    "${CONDA}" env create -p "${prefix}" -f "${REPO_DIR}/envs/$1.yaml"
  md5sum < "${REPO_DIR}/envs/$1.yaml" > "${prefix}/.schla_ready"
}

check_snakemake() {
  env_ready snakemake \
    && "${RES}/envs/snakemake/bin/snakemake" --version >/dev/null 2>&1
}

install_snakemake() { install_env snakemake; }

check_arcashla() {
  local bin="${RES}/envs/arcashla/bin"
  env_ready arcashla \
    && [ "$(cat "${RES}/envs/arcashla/.schla_hla_reference" 2>/dev/null)" = "${IMGT_COMMIT}" ] \
    && [ -x "${bin}/arcasHLA" ] && "${bin}/samtools" --version >/dev/null 2>&1 \
    && grep -qw pseudo <<< "$("${bin}/kallisto" 2>&1)"
}

install_arcashla() {
  local prefix="${RES}/envs/arcashla" bin="${RES}/envs/arcashla/bin" imgt
  install_env arcashla
  imgt=$(newest "${prefix}"/share/arcas-hla-*/dat)/IMGTHLA
  log "downloading the IMGT/HLA database (release ${IMGT_RELEASE}) for arcasHLA"
  # arcasHLA clones the IMGT/HLA repository (its own build of the latest release can
  # fail, see below), then the tested release is checked out. IMGT/HLA now ships
  # hla.dat zipped, which arcasHLA 0.6.0 does not expect: unzip it, then build
  # arcasHLA's reference from it.
  PATH="${bin}:${PATH}" "${bin}/arcasHLA" reference --update || true
  if [ ! -d "${imgt}/.git" ]; then
    echo "The IMGT/HLA download failed: no git clone in ${imgt}" >&2
    return 1
  fi
  rm -f "${imgt}/hla.dat"
  PATH="${bin}:${PATH}" git -C "${imgt}" checkout -q -f "${IMGT_COMMIT}"
  if [ ! -f "${imgt}/hla.dat" ]; then
    if [ ! -f "${imgt}/hla.dat.zip" ]; then
      echo "The IMGT/HLA release ${IMGT_RELEASE} has no hla.dat or hla.dat.zip in ${imgt}" >&2
      return 1
    fi
    "${bin}/python" -c 'import sys, zipfile; zipfile.ZipFile(sys.argv[1]).extract("hla.dat", sys.argv[2])' \
      "${imgt}/hla.dat.zip" "${imgt}"
  fi
  PATH="${bin}:${PATH}" "${bin}/arcasHLA" reference --rebuild
  echo "${IMGT_COMMIT}" > "${prefix}/.schla_hla_reference"
}

# --------------------------------------------------------------------------

STEPS=(
  "cellranger|Cell Ranger"
  "vdj|10x V(D)J reference (${VDJ_REF})"
  "gex|10x GEX reference (${GEX_REF})"
  "subset_bam|subset-bam"
  "snakemake|Snakemake"
  "arcashla|arcasHLA + samtools + IMGT/HLA database"
  "apptainer|apptainer (runs the souporcell image)"
  "souporcell|souporcell image"
)

mkdir -p "${RES}"

# What the machine must already have; setup.sh installs everything else.
if ! ${CHECK_ONLY}; then
  missing=()
  for tool in curl git tar; do
    command -v "${tool}" >/dev/null || missing+=("${tool}")
  done
  if [ ${#missing[@]} -gt 0 ]; then
    echo "Please install ${missing[*]} first (or ask your system administrator), then run ./setup.sh again." >&2
    exit 1
  fi
fi
failed=()
summary=()
for step in "${STEPS[@]}"; do
  key=${step%%|*}
  label=${step#*|}
  if "check_${key}"; then
    summary+=("  ok       ${label}")
    continue
  fi
  if ${CHECK_ONLY}; then
    summary+=("  MISSING  ${label}")
    failed+=("${label}")
    continue
  fi
  log "=== ${label}"
  set +e
  ( set -e; "install_${key}" )
  status=$?
  set -e
  if [ ${status} -eq 0 ] && "check_${key}"; then
    summary+=("  ok       ${label} (installed now)")
  else
    summary+=("  FAILED   ${label}")
    failed+=("${label}")
  fi
done

echo
echo "scHLA-demux setup (${RES}):"
printf '%s\n' "${summary[@]}"
echo
if [ ${#failed[@]} -gt 0 ]; then
  if ${CHECK_ONLY}; then
    echo "Run ./setup.sh to install what is missing."
  else
    echo "Some parts are not installed: see the messages above, fix them and run ./setup.sh again."
  fi
  exit 1
fi
echo "Everything is installed. Next: fill in config/samples.csv and run the dry run (see README)."
