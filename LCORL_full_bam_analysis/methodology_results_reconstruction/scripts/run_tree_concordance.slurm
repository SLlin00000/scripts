#!/bin/bash
#SBATCH --job-name=LCORL_tree_distance
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/LCORL_tree_distance_%j.out
#SBATCH --error=logs/LCORL_tree_distance_%j.err

set -euo pipefail

WORK=/QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis/methodology_results_reconstruction
PROJECT=/QRISdata/Q9774/lcorl_project

cd "$WORK"

RESULT="$WORK/results/tree_concordance_${SLURM_JOB_ID}"
mkdir "$RESULT"

python scripts/analyse_tree_concordance.py \
  --project "$PROJECT" \
  --out "$RESULT" \
  --permutations 9999 \
  --seed 7621

echo "PASS" > "$RESULT/COMPLETED_VALIDATED.txt"

echo "===== ANALYSIS COMPLETED ====="
echo "Results: $RESULT"
