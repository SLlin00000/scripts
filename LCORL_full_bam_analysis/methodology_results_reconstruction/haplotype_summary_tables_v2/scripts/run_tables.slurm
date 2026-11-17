#!/bin/bash
#SBATCH --job-name=LCORL_tables_v2
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=presentation_tables_v2/logs/tables_%j.out
#SBATCH --error=presentation_tables_v2/logs/tables_%j.err

set -euo pipefail

cd /QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis/methodology_results_reconstruction

python presentation_tables_v2/scripts/build_tables.py
