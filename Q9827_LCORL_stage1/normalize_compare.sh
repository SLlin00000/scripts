#!/bin/bash
#SBATCH --job-name=LCORL_compare
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/user/s4922196/Q9827_LCORL_comparison/compare.out
#SBATCH --error=/scratch/user/s4922196/Q9827_LCORL_comparison/compare.err

module load bcftools/1.18-gcc-12.3.0

REF=/scratch/user/s4922196/Q9827_reference/ARS-UCD1.2_Btau5.0.1Y.fa
Q9774=/scratch/user/s4922196/Q9827_LCORL_comparison/Q9774_Run9_LCORL.vcf.gz
Q9827=/scratch/user/s4922196/Q9827_LCORL_joint/Q9827_LCORL_joint.vcf.gz

bcftools norm \
  -f "$REF" \
  -m -any \
  "$Q9774" \
  -Ou \
| bcftools view -v snps -Oz \
  -o Q9774_Run9_LCORL_snps.norm.vcf.gz

bcftools index -t Q9774_Run9_LCORL_snps.norm.vcf.gz

bcftools norm \
  -f "$REF" \
  -m -any \
  "$Q9827" \
  -Ou \
| bcftools view -v snps -Oz \
  -o Q9827_LCORL_snps.norm.vcf.gz

bcftools index -t Q9827_LCORL_snps.norm.vcf.gz

mkdir -p isec

bcftools isec \
  -p isec \
  Q9774_Run9_LCORL_snps.norm.vcf.gz \
  Q9827_LCORL_snps.norm.vcf.gz

echo "Q9774 normalized SNPs: $(bcftools view -H Q9774_Run9_LCORL_snps.norm.vcf.gz | wc -l)"
echo "Q9827 normalized SNPs: $(bcftools view -H Q9827_LCORL_snps.norm.vcf.gz | wc -l)"

for f in isec/000*.vcf; do
    echo "$(basename "$f"): $(grep -vc '^#' "$f")"
done

echo "Finished successfully"

