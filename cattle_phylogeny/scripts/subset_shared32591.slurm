#!/bin/bash
#SBATCH --job-name=shared32591
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/shared32591_%j.out
#SBATCH --error=logs/shared32591_%j.err

set -euo pipefail

BASE="/QRISdata/Q9774/lcorl_project/cattle_phylogeny"

module load bcftools/1.22-gcc-14.2.0

REGIONS="$BASE/metadata/shared_exact_clean_positions.tsv"

LOCAL_IN="$BASE/results/local174_PASS_biallelicSNPs_LCORL.vcf.gz"
LOCAL_OUT="$BASE/results/local174_shared32591.vcf.gz"

RUN8_IN="$BASE/run8/Run8_1089_KEEP_biallelicSNPs_LCORL.vcf.gz"
RUN8_OUT="$BASE/results/run8_1089_shared32591.vcf.gz"

echo "Start: $(date)"

echo "Subsetting LOCAL..."
bcftools view \
    -R "$REGIONS" \
    -Oz \
    -o "$LOCAL_OUT" \
    "$LOCAL_IN"

bcftools index -f -t "$LOCAL_OUT"

echo "Subsetting RUN8..."
bcftools view \
    -R "$REGIONS" \
    -Oz \
    -o "$RUN8_OUT" \
    "$RUN8_IN"

bcftools index -f -t "$RUN8_OUT"

echo "LOCAL samples:"
bcftools query -l "$LOCAL_OUT" | wc -l

echo "LOCAL variants:"
bcftools view -H "$LOCAL_OUT" | wc -l

echo "RUN8 samples:"
bcftools query -l "$RUN8_OUT" | wc -l

echo "RUN8 variants:"
bcftools view -H "$RUN8_OUT" | wc -l

echo "Finished: $(date)"

