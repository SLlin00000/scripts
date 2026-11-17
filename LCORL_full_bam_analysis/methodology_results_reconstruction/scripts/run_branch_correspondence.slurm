#!/bin/bash
#SBATCH --job-name=LCORL_branch_map
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/LCORL_branch_map_%j.out
#SBATCH --error=logs/LCORL_branch_map_%j.err

set -euo pipefail

WORK=/QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis/methodology_results_reconstruction
PROJECT=/QRISdata/Q9774/lcorl_project

cd "$WORK"

RESULT="$WORK/results/branch_mapping_${SLURM_JOB_ID}"
mkdir "$RESULT"

python scripts/build_branch_correspondence.py \
  --project "$PROJECT" \
  --out "$RESULT"

echo PASS > "$RESULT/COMPLETED_VALIDATED.txt"

echo "===== BRANCH MAPPING COMPLETED ====="
echo "$RESULT"
