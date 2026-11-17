#!/bin/bash -l

#SBATCH --job-name=VEP107_cache
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/VEP107_cache_%j.out
#SBATCH --error=logs/VEP107_cache_%j.err

set -euo pipefail

PROJECT_ROOT="/QRISdata/Q9774/lcorl_project"
DOWNLOAD_DIR="${PROJECT_ROOT}/resources/downloads"
CACHE_DIR="${PROJECT_ROOT}/resources/vep_cache"

CACHE_URL="https://ftp.ensembl.org/pub/release-107/variation/indexed_vep_cache/bos_taurus_vep_107_ARS-UCD1.2.tar.gz"
TARBALL="${DOWNLOAD_DIR}/bos_taurus_vep_107_ARS-UCD1.2.tar.gz"

mkdir -p "$DOWNLOAD_DIR" "$CACHE_DIR"

echo "Node: $(hostname)"
echo "Start: $(date)"
echo "Downloading to: $TARBALL"

wget \
  --continue \
  --tries=5 \
  --timeout=60 \
  --output-document="$TARBALL" \
  "$CACHE_URL"

echo "Checking downloaded archive..."
gzip -t "$TARBALL"

echo "Extracting cache..."
tar -xzf "$TARBALL" -C "$CACHE_DIR"

echo "Installed cache directories:"
find "$CACHE_DIR/bos_taurus" -maxdepth 2 -type d -print

echo "Cache size:"
du -sh "$CACHE_DIR/bos_taurus"

echo "Finish: $(date)"

