#!/bin/bash
#SBATCH --job-name=run8_KEEP
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/run8_KEEP_%j.out
#SBATCH --error=logs/run8_KEEP_%j.err

set -euo pipefail

BASE="/QRISdata/Q9774/lcorl_project/cattle_phylogeny"

module load bcftools/1.22-gcc-14.2.0

IN="$BASE/run8/Run8_Chr6_LCORL_plusminus1Mb.vcf.gz"
SAMPLES="$BASE/metadata/run8_KEEP_clean_samples.txt"
OUT="$BASE/run8/Run8_1089_KEEP_biallelicSNPs_LCORL.vcf.gz"

echo "Start: $(date)"

bcftools view \
    -S "$SAMPLES" \
    -m2 -M2 \
    -v snps \
    -Oz \
    -o "$OUT" \
    "$IN"

bcftools index -f -t "$OUT"

echo "Samples:"
bcftools query -l "$OUT" | wc -l

echo "Variants:"
bcftools view -H "$OUT" | wc -l

echo "Output:"
ls -lh "$OUT" "$OUT.tbi"

echo "Finished: $(date)"

