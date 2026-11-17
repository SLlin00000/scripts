#!/bin/bash
#SBATCH --job-name=LCORL_blocks
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis/haplotype/blocks_%j.out
#SBATCH --error=/QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis/haplotype/blocks_%j.err

set -euo pipefail

cd /QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis

python3 haplotype/find_candidate_blocks.py

echo "Candidate block discovery finished successfully"

