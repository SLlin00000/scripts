#!/bin/bash

MAP="/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/run8_sample_accession_map.tsv"
OUT="/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/test10_metadata.tsv"

echo -e "ENA_accession\tBioSample\tisolate\torganism\tsex\tdescription" > "$OUT"

head -10 "$MAP" | while IFS=$'\t' read -r ena biosample
do
    tmp=$(mktemp)

    curl -s \
    "https://www.ebi.ac.uk/ena/browser/api/xml/$ena" \
    > "$tmp"

    isolate=$(awk '
        /<TAG>isolate<\/TAG>/ {getline; gsub(/.*<VALUE>|<\/VALUE>.*/, ""); print}
    ' "$tmp")

    organism=$(awk '
        /<TAG>organism<\/TAG>/ {getline; gsub(/.*<VALUE>|<\/VALUE>.*/, ""); print}
    ' "$tmp")

    sex=$(awk '
        /<TAG>sex<\/TAG>/ {getline; gsub(/.*<VALUE>|<\/VALUE>.*/, ""); print}
    ' "$tmp")

    description=$(sed -n 's:.*<DESCRIPTION>\(.*\)</DESCRIPTION>.*:\1:p' "$tmp")

    echo -e "${ena}\t${biosample}\t${isolate}\t${organism}\t${sex}\t${description}" >> "$OUT"

    rm -f "$tmp"
    sleep 1
done

