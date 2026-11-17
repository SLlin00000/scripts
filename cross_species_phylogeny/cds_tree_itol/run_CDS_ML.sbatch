#!/bin/bash
#SBATCH --job-name=LCORL_CDS_ML
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/QRISdata/Q9774/lcorl_project/cross_species_phylogeny/cds_tree_itol/CDS_ML-%j.out

set -euo pipefail

module load mafft/7.526-gcc-13.3.0-with-extensions
module load iqtree/2.2.2.3--h2202e69_2/module
export SINGULARITY_BINDPATH="${SINGULARITY_BINDPATH:+${SINGULARITY_BINDPATH},}/QRISdata"

cd /QRISdata/Q9774/lcorl_project/cross_species_phylogeny

mafft --nuc --auto --thread "${SLURM_CPUS_PER_TASK}" \
  cds_tree_itol/LCORL_CDS_12taxa.fna \
  > cds_tree_itol/LCORL_CDS_12taxa.mafft.fasta \
  2> cds_tree_itol/LCORL_CDS_12taxa.mafft.log

test -s cds_tree_itol/LCORL_CDS_12taxa.mafft.fasta

iqtree2 \
  -s cds_tree_itol/LCORL_CDS_12taxa.mafft.fasta \
  -st DNA \
  -m MFP \
  -b 1000 \
  -nt "${SLURM_CPUS_PER_TASK}" \
  -seed 20261021 \
  -o horse \
  -pre cds_tree_itol/LCORL_CDS_12taxa_ML_standardBS1000
