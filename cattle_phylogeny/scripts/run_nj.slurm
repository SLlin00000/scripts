#!/bin/bash
#SBATCH --job-name=LCORL_NJ
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/LCORL_NJ_%j.out
#SBATCH --error=logs/LCORL_NJ_%j.err

cd /QRISdata/Q9774/lcorl_project/cattle_phylogeny

./scripts/neighbor_joining

