#!/bin/bash
#SBATCH --job-name=LCORL_IsoSeq_ovary
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

module load minimap2/2.26-gcccore-12.3.0
module load samtools/1.18-gcc-12.3.0

REF=/scratch/user/s4922196/Q9827_reference/ARS-UCD1.2_Btau5.0.1Y.fa
READS=/QRISdata/Q9774/lcorl_project/isoseq/Ovary_minpass3_nopolish_ccs_lima.5p--3p.flnc.fastq

test -s "$REF"
test -s "$READS"

minimap2 --version
samtools --version | head -n 1

# 比对整个参考基因组；暂不强制转录本方向。
minimap2 -t 6 -ax splice:hq "$REF" "$READS" \
  2> minimap2.log \
  | samtools sort -@ 1 -m 2G -o Ovary.genome.sorted.bam -

samtools quickcheck -v Ovary.genome.sorted.bam
samtools index Ovary.genome.sorted.bam
samtools flagstat Ovary.genome.sorted.bam > Ovary.flagstat.txt

# 确认参考中的 chr6 名称。
CHR=$(samtools idxstats Ovary.genome.sorted.bam \
  | awk '$1=="6" || $1=="chr6" {print $1}')

if [[ -z "$CHR" || "$CHR" == *$'\n'* ]]; then
    printf 'ERROR: Cannot uniquely identify chromosome 6.\n' >&2
    exit 1
fi

REGION="${CHR}:37377806-37556970"

# 提取与 LCORL 区域重叠的主要比对，暂不按 MAPQ 筛除。
samtools view -b -F 2308 \
  Ovary.genome.sorted.bam "$REGION" \
  > Ovary.LCORL.candidates.bam

samtools index Ovary.LCORL.candidates.bam

samtools fasta -n Ovary.LCORL.candidates.bam \
  > Ovary.LCORL.candidates.fasta

samtools view Ovary.LCORL.candidates.bam \
  | awk 'BEGIN {
      OFS="\t"; print "read_id","flag","chrom","position","MAPQ","CIGAR"
    }
    {print $1,$2,$3,$4,$5,$6}' \
  > Ovary.LCORL.alignments.tsv

{
    printf 'Region: %s\n' "$REGION"
    printf 'LCORL-overlapping primary alignments: '
    samtools view -c Ovary.LCORL.candidates.bam
    printf 'Candidate FASTA records: '
    awk '/^>/ {n++} END {print n+0}' Ovary.LCORL.candidates.fasta
} > summary.txt

cat summary.txt
printf 'DONE\n'
