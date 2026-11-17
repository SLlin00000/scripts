#!/bin/bash
#SBATCH --job-name=Q9827_LCORL
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --array=1-125%10
#SBATCH --output=/scratch/user/s4922196/Q9827_LCORL/logs/%A_%a.out
#SBATCH --error=/scratch/user/s4922196/Q9827_LCORL/logs/%A_%a.err

module load samtools/1.18-gcc-12.3.0

LIST=/scratch/user/s4922196/Q9827_LCORL/bam_list.txt
OUTDIR=/scratch/user/s4922196/Q9827_LCORL/regional_bam
REGION="6:36377806-38557145"

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" )

if [ -z "$BAM" ]; then
	 echo "ERROR: No BAM found for array task ${SLURM_ARRAY_TASK_ID}" 
	exit 1
fi
BASE=$(basename "$BAM" .bam)
OUT="${OUTDIR}/${BASE}.LCORL_plusminus1Mb.bam"

echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "INPUT: $BAM"
echo "Region: $REGION"
echo "Output: $OUT"

samtools quickcheck -v "$BAM"
if [ $? -ne 0 ]; then
	echo "ERROR: input BAM faild quickcheck"
	exit 1
fi

samtools view\
	-b \
	"$BAM" \
	"$REGION" \
	-o "$OUT"

if [ $? -ne 0 ]; then
	echo "ERROR: samtools view failed"
	exit 1
fi

samtools index "$OUT" 

if [ $? -ne 0 ]; then
	echo "ERROR: samtools index failed"
	exit 1
fi

samtools quickcheck -v "$OUT"
if [ $? -ne 0 ]; then
	echo "ERROR: output BAM failed quickcheck"
	eexit 1
fi

COUNT=$(samtools view -c "$OUT")

echo "Alignment records: $COUNT"
echo "Finished successfully"


