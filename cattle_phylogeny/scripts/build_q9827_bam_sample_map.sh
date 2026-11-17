#!/bin/bash

out="metadata/q9827_bam_to_sample.tsv"

echo -e "bam_file\tsample_id" > "$out"

find /QRISdata/Q9827/BAMandBAI -type f -name '*.bam' | sort | while read -r bam
do
    sm=$(samtools view -H "$bam" \
        | awk -F'\t' '$1=="@RG" {
            for(i=1;i<=NF;i++){
                if($i ~ /^SM:/){
                    sub(/^SM:/,"",$i)
                    print $i
                    exit
                }
            }
        }')

    echo -e "$(basename "$bam")\t$sm" >> "$out"
done
