#!/bin/bash
#SBATCH --job-name=LCORL_correspondence
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/LCORL_correspondence_%j.out
#SBATCH --error=logs/LCORL_correspondence_%j.err

set -euo pipefail

WORK=/QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis/methodology_results_reconstruction
PROJECT=/QRISdata/Q9774/lcorl_project

cd "${WORK}"

JOB_ID="${SLURM_JOB_ID}"
RESULT="${WORK}/results/run_${JOB_ID}"

mkdir "${RESULT}"

echo "===== LCORL RECONSTRUCTION ====="
echo "Job ID: ${JOB_ID}"
echo "Started: $(date)"
echo "Output: ${RESULT}"

python "${WORK}/scripts/reconstruct_correspondence.py" \
    --project "${PROJECT}" \
    --out "${RESULT}"

echo "PASS" > "${RESULT}/COMPLETED_VALIDATED.txt"

echo "===== JOB COMPLETED ====="
echo "Finished: $(date)"
echo "Results: ${RESULT}"
