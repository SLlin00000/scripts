#!/bin/bash
#SBATCH --job-name=LCORL_AA_ML
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/QRISdata/Q9774/lcorl_project/cross_species_phylogeny/aa_tree_itol/AA_ML-%j.out

set -euo pipefail

module load mafft/7.526-gcc-13.3.0-with-extensions
module load iqtree/2.2.2.3--h2202e69_2/module
export SINGULARITY_BINDPATH="${SINGULARITY_BINDPATH:+${SINGULARITY_BINDPATH},}/QRISdata"

cd /QRISdata/Q9774/lcorl_project/cross_species_phylogeny

mafft --amino --auto --thread "${SLURM_CPUS_PER_TASK}" \
  aa_tree_itol/LCORL_AA_12taxa.faa \
  > aa_tree_itol/LCORL_AA_12taxa.mafft.faa \
  2> aa_tree_itol/LCORL_AA_12taxa.mafft.log

test -s aa_tree_itol/LCORL_AA_12taxa.mafft.faa

iqtree2 \
  -s aa_tree_itol/LCORL_AA_12taxa.mafft.faa \
  -st AA \
  -m MFP \
  -b 1000 \
  -nt "${SLURM_CPUS_PER_TASK}" \
  -seed 20261021 \
  -o horse \
  -pre aa_tree_itol/LCORL_AA_12taxa_ML_standardBS1000
