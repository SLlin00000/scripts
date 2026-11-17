#!/bin/bash
#SBATCH --job-name=run8_meta10
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/run8_meta10_%j.out
#SBATCH --error=logs/run8_meta10_%j.err

set -euo pipefail

MAP="/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/run8_sample_accession_map.tsv"
OUT="/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/test10_metadata.tsv"

echo -e "ENA_accession\tBioSample\tisolate\torganism\tsex\tdescription" > "$OUT"

head -10 "$MAP" | while IFS=$'\t' read -r ena biosample
do
    tmp=$(mktemp)

    echo "Processing $ena $biosample"

    curl --connect-timeout 20 \
         --max-time 60 \
         --retry 3 \
         --retry-delay 3 \
         -s \
         "https://www.ebi.ac.uk/ena/browser/api/xml/$ena" \
         > "$tmp"

    isolate=$(awk '
        /<TAG>isolate<\/TAG>/ {
            getline
            gsub(/.*<VALUE>|<\/VALUE>.*/, "")
            print
            exit
        }
    ' "$tmp")

    organism=$(awk '
        /<TAG>organism<\/TAG>/ {
            getline
            gsub(/.*<VALUE>|<\/VALUE>.*/, "")
            print
            exit
        }
    ' "$tmp")

    sex=$(awk '
        /<TAG>sex<\/TAG>/ {
            getline
            gsub(/.*<VALUE>|<\/VALUE>.*/, "")
            print
            exit
        }
    ' "$tmp")

    description=$(sed -n \
        's:.*<DESCRIPTION>\(.*\)</DESCRIPTION>.*:\1:p' \
        "$tmp" | head -1)

    echo -e "${ena}\t${biosample}\t${isolate}\t${organism}\t${sex}\t${description}" >> "$OUT"

    rm -f "$tmp"
    sleep 1
done

echo "Finished"
column -t -s $'\t' "$OUT"

