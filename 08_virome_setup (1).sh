#!/usr/bin/env bash
#==============================================================================
# 08_virome_setup.sh
# One-time installation of the software environment and Kraken2 viral database
# for 08_virome_pipeline.sh (Ubuntu / WSL2)
#
# Authors : <AUTHOR NAMES>
# Paper   : <CITATION / DOI OF PUBLISHED ARTICLE>
# Licence : <LICENCE, e.g. MIT>
#
# Installs Miniconda (if missing), a conda environment with fastp, Kraken2,
# seqkit, pigz, R and the R packages used by 08_virome_analysis.R, and the
# Kraken2 RefSeq viral database (one fixed release, see KRAKEN_DB_URL).
#
# For exact reproducibility, export the environment after installation and
# deposit the file with the scripts:
#   conda env export -n potato_virome > environment.yml
#==============================================================================
set -euo pipefail

# -------------------- SETTINGS --------------------
ENV_NAME="potato_virome"
DB_DIR="<PATH_TO_KRAKEN2_VIRAL_DB_FOLDER>"       # e.g. $HOME/kraken2_dbs/k2_viral
# One fixed database release, so that every installation classifies against
# the same reference (see https://benlangmead.github.io/aws-indexes/k2)
KRAKEN_DB_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240904.tar.gz"

if [[ "$DB_DIR" == *"<"* ]]; then
  echo "ERROR: set DB_DIR at the top of this script before running." >&2; exit 1
fi

echo "============================================================"
echo " [1/4] Checking system utilities and Miniconda"
echo "============================================================"
command -v curl >/dev/null 2>&1 || { echo "Installing curl..."; sudo apt-get update && sudo apt-get install -y curl; }
command -v tar  >/dev/null 2>&1 || { echo "Installing tar...";  sudo apt-get update && sudo apt-get install -y tar; }

if [ ! -d "$HOME/miniconda3" ]; then
  echo "Installing Miniconda to ~/miniconda3..."
  tmp=$(mktemp -d)
  curl -fsSL -o "$tmp/Miniconda3.sh" https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash "$tmp/Miniconda3.sh" -b -p "$HOME/miniconda3"
  rm -rf "$tmp"
else
  echo "Miniconda already exists at ~/miniconda3 — skipping installation."
fi
# shellcheck disable=SC1091
source "$HOME/miniconda3/etc/profile.d/conda.sh"

echo ""
echo "============================================================"
echo " [2/4] Conda environment: ${ENV_NAME}"
echo "============================================================"
if conda info --envs | grep -E "^${ENV_NAME}\s" >/dev/null 2>&1; then
  echo "Environment '${ENV_NAME}' already exists — skipping creation."
else
  conda create -y -n "${ENV_NAME}" python=3.10
fi
conda activate "${ENV_NAME}"

echo ""
echo "============================================================"
echo " [3/4] Installing required packages"
echo "============================================================"
CONDA_PKGS=(fastp kraken2 seqkit pigz r-base r-vegan r-ggplot2 r-dplyr r-tidyr
            r-readr r-tibble r-scales r-pheatmap r-ape)
MISSING_PKGS=()
for pkg in "${CONDA_PKGS[@]}"; do
  conda list | grep -E "^${pkg}\s" >/dev/null 2>&1 || MISSING_PKGS+=("${pkg}")
done
if [ ${#MISSING_PKGS[@]} -eq 0 ]; then
  echo "All packages are already installed."
else
  echo "Installing: ${MISSING_PKGS[*]}"
  conda install -y -c conda-forge -c bioconda "${MISSING_PKGS[@]}"
fi

echo ""
echo "Installed versions:"
fastp --version 2>&1 | head -1 || true
kraken2 --version 2>&1 | head -1 || true
Rscript -e 'cat("R:", R.version.string, "\n")'
Rscript -e 'suppressPackageStartupMessages(library(vegan)); cat("vegan:", as.character(packageVersion("vegan")), "\n")'

echo ""
echo "============================================================"
echo " [4/4] Kraken2 viral database"
echo "============================================================"
if [ -f "${DB_DIR}/hash.k2d" ]; then
  echo "Kraken2 database present in: ${DB_DIR}"
else
  echo "Downloading $(basename "${KRAKEN_DB_URL}")..."
  mkdir -p "${DB_DIR}"
  tmp=$(mktemp -d)
  curl -fsSL -o "$tmp/k2_viral.tar.gz" "${KRAKEN_DB_URL}"
  tar -xzf "$tmp/k2_viral.tar.gz" -C "${DB_DIR}"
  find "${DB_DIR}" -mindepth 2 -type f -exec mv -t "${DB_DIR}" {} +
  rm -rf "$tmp"
  basename "${KRAKEN_DB_URL}" > "${DB_DIR}/DATABASE_VERSION.txt"
  echo "Database unpacked to ${DB_DIR}"
fi

echo ""
echo "============================================================"
echo " Setup complete. Next: bash 08_virome_pipeline.sh"
echo "============================================================"
