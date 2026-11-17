#!/bin/bash
#SBATCH --job-name=LCORL_gene_ML
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/QRISdata/Q9774/lcorl_project/cross_species_phylogeny/gene_tree_itol/gene_ML-%j.out

set -euo pipefail

module load iqtree/2.2.2.3--h2202e69_2/module

export SINGULARITY_BINDPATH="${SINGULARITY_BINDPATH:+${SINGULARITY_BINDPATH},}/QRISdata"

cd /QRISdata/Q9774/lcorl_project/cross_species_phylogeny

iqtree2 \
  -s alignments/LCORL_gene_12taxa.mafft.fasta \
  -st DNA \
  -m MFP \
  -b 1000 \
  -nt "${SLURM_CPUS_PER_TASK}" \
  -seed 20261020 \
  -o horse \
  -pre gene_tree_itol/LCORL_gene_12taxa_ML_standardBS1000
