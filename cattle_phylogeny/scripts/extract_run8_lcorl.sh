#!/bin/bash
#SBATCH --job-name=run8_lcorl
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/run8_lcorl_%j.out
#SBATCH --error=logs/run8_lcorl_%j.err

module load bcftools/1.22-gcc-14.2.0

set -euo pipefail

URL="https://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ173/ERZ1738264/Chr6-Run8-TAUIND-public.vcf.gz"

OUT="/QRISdata/Q9774/lcorl_project/cattle_phylogeny/run8/Run8_Chr6_LCORL_plusminus1Mb.vcf.gz"

mkdir -p /QRISdata/Q9774/lcorl_project/cattle_phylogeny/run8

echo "Start: $(date)"

bcftools view \
-r 6:36377806-38557145 \
-Oz \
-o "$OUT" \
"$URL"

bcftools index -t "$OUT"

echo "Samples:"
bcftools query -l "$OUT" | wc -l

echo "Variants:"
bcftools index -n "$OUT"

ls -lh "$OUT" "$OUT.tbi"

echo "Finished: $(date)"
