#!/bin/bash
#SBATCH --job-name=run8_meta
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --array=0-9%5
#SBATCH --output=logs/run8_meta_%A_%a.out
#SBATCH --error=logs/run8_meta_%A_%a.err

set -u

BASE="/QRISdata/Q9774/lcorl_project/cattle_phylogeny"
MAP="$BASE/metadata/run8_sample_accession_map.tsv"
PARSER="$BASE/scripts/parse_ena_sample.py"
PARTDIR="$BASE/metadata/run8_parts"
XMLDIR="$BASE/metadata/run8_xml"

mkdir -p "$PARTDIR" "$XMLDIR"

TOTAL=$(wc -l < "$MAP")
NTASKS=10
CHUNK=$(( (TOTAL + NTASKS - 1) / NTASKS ))

START=$(( SLURM_ARRAY_TASK_ID * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

if [ "$END" -gt "$TOTAL" ]; then
    END=$TOTAL
fi

OUT="$PARTDIR/run8_metadata_part_${SLURM_ARRAY_TASK_ID}.tsv"

echo -e "ENA_accession\tBioSample\tscientific_name\torganism_attribute\tisolate\tbreed\tsex\tcountry\tsample_name\tdev_stage\tdescription" > "$OUT"

echo "Task $SLURM_ARRAY_TASK_ID"
echo "Rows $START-$END"

sed -n "${START},${END}p" "$MAP" |
while IFS=$'\t' read -r ena biosample
do
    [ -z "$ena" ] && continue

    xml="$XMLDIR/${ena}.xml"

    echo "Fetching $ena $biosample"

    if ! curl \
        --connect-timeout 20 \
        --max-time 60 \
        --retry 3 \
        --retry-delay 2 \
        -fsS \
        "https://www.ebi.ac.uk/ena/browser/api/xml/$ena" \
        -o "$xml"
    then
        echo -e "${ena}\t${biosample}\tFETCH_FAILED\t\t\t\t\t\t\t\t" >> "$OUT"
        rm -f "$xml"
        continue
    fi

    parsed=$(python3 "$PARSER" "$xml")

    echo -e "${ena}\t${biosample}\t${parsed}" >> "$OUT"

    sleep 0.5
done

echo "Finished task $SLURM_ARRAY_TASK_ID"
wc -l "$OUT"

