#!/bin/bash
#SBATCH --job-name=LCORL_IsoSeq_screen
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

module load minimap2/2.26-gcccore-12.3.0
module load samtools/1.18-gcc-12.3.0

ROOT="$PWD"
LINE=$((SLURM_ARRAY_TASK_ID + 1))
READS=$(sed -n "${LINE}p" "$ROOT/samples.txt")
test -s "$READS"

SAMPLE=$(basename "$READS" .flnc.fastq)
mkdir -p "$ROOT/$SAMPLE"
cd "$ROOT/$SAMPLE"

REF=/scratch/user/s4922196/Q9827_reference/ARS-UCD1.2_Btau5.0.1Y.fa
test -s "$REF"

minimap2 -t 6 -ax splice:hq "$REF" "$READS" \
  2> minimap2.log \
  | samtools sort -@ 1 -m 2G -o genome.sorted.bam -

samtools quickcheck -v genome.sorted.bam
samtools index genome.sorted.bam
samtools flagstat genome.sorted.bam > flagstat.txt

CHR=$(samtools idxstats genome.sorted.bam \
  | awk '$1=="6" || $1=="chr6" {print $1}')

if [[ -z "$CHR" || "$CHR" == *$'\n'* ]]; then
    printf 'ERROR: Cannot uniquely identify chromosome 6.\n' >&2
    exit 1
fi

GENE="${CHR}:37380296-37557106"
FLANK="${CHR}:37370296-37567106"

# 保存准确基因区间内所有已比对记录，供后续核验。
samtools view -b -F 4 genome.sorted.bam "$GENE" \
  > LCORL.all_alignments.bam
samtools index LCORL.all_alignments.bam

# 主要比对对应的候选序列，尚不是最终 isoform。
samtools view -b -F 2308 genome.sorted.bam "$GENE" \
  > LCORL.primary.bam
samtools index LCORL.primary.bam
samtools fasta -n LCORL.primary.bam > LCORL.primary.fasta

PRIMARY=$(samtools view -c -F 2308 genome.sorted.bam "$GENE")
SECONDARY=$(samtools view -c -f 256 -F 4 genome.sorted.bam "$GENE")
SUPPLEMENTARY=$(samtools view -c -f 2048 -F 4 genome.sorted.bam "$GENE")
FLANK_ALL=$(samtools view -c -F 4 genome.sorted.bam "$FLANK")

{
    printf 'sample\tprimary\tsecondary\tsupplementary\tflank_all\n'
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$SAMPLE" "$PRIMARY" "$SECONDARY" "$SUPPLEMENTARY" "$FLANK_ALL"
} > summary.tsv

cat summary.tsv
printf 'DONE\n'
