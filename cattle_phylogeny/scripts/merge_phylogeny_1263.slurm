#!/bin/bash
#SBATCH --job-name=merge1263
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/merge1263_%j.out
#SBATCH --error=logs/merge1263_%j.err

set -euo pipefail

BASE="/QRISdata/Q9774/lcorl_project/cattle_phylogeny"

module load bcftools/1.22-gcc-14.2.0

LOCAL="$BASE/results/local174_shared32591.vcf.gz"
RUN8="$BASE/results/run8_1089_shared32591.vcf.gz"
OUT="$BASE/results/cattle1263_shared32591_merged.vcf.gz"

echo "Start: $(date)"

bcftools merge \
    "$LOCAL" \
    "$RUN8" \
    -Oz \
    -o "$OUT"

bcftools index -f -t "$OUT"

echo "Samples:"
bcftools query -l "$OUT" | wc -l

echo "Variants:"
bcftools view -H "$OUT" | wc -l

echo "Output:"
ls -lh "$OUT" "$OUT.tbi"

echo "Finished: $(date)"

