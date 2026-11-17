#!/bin/bash
#SBATCH --job-name=LCORL_phase
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis/haplotype/beagle_%j.out
#SBATCH --error=/QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis/haplotype/beagle_%j.err

set -euo pipefail

module purge
module load beagle/5.4.22jul22.46e-java-11

WORK=/QRISdata/Q9774/lcorl_project/LCORL_full_bam_analysis

echo "Beagle jar: ${EBROOTBEAGLE}/beagle.jar"
echo "Input: ${WORK}/haplotype/LCORL_QC_for_phasing.vcf.gz"

java -Xmx12g -jar ${EBROOTBEAGLE}/beagle.jar \
  gt=${WORK}/haplotype/LCORL_QC_for_phasing.vcf.gz \
  out=${WORK}/haplotype/LCORL_beagle_phased \
  nthreads=4

echo "Beagle phasing finished successfully"

