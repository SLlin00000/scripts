#!/bin/bash -l

#SBATCH --job-name=LCORL_interpro
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/LCORL_interpro_%j.out
#SBATCH --error=logs/LCORL_interpro_%j.err

set -euo pipefail

cd /QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis

mkdir -p protein/interproscan logs

module purge
module load interproscan/5.76-107.0

echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "InterProScan command: $(command -v interproscan)"

interproscan \
    -i protein/LCORL_ENSBTAP00000060198.fa \
    -b protein/interproscan/LCORL_interpro \
    -f TSV,GFF3,JSON \
    -t p \
    -iprlookup \
    -goterms \
    -pa \
    -dp \
    -cpu 8

echo "Finish time: $(date)"
