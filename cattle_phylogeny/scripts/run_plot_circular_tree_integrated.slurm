#!/bin/bash
#SBATCH --job-name=tree_plot
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/tree_plot_%j.out
#SBATCH --error=logs/tree_plot_%j.err

cd /QRISdata/Q9774/lcorl_project/cattle_phylogeny

python3 scripts/plot_circular_tree_integrated.py
