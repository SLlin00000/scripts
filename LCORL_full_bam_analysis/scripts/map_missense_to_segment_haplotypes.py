#!/usr/bin/env python3

import csv
import gzip
from collections import defaultdict
from pathlib import Path

vcf_file = Path("haplotype/LCORL_gene_phased_MAF05.vcf.gz")
segment_file = Path("haplotype/segment_SNP_annotation.tsv")
evidence_file = Path(
    "summary/functional_annotation/"
    "LCORL_missense_integrated_evidence.tsv"
)
output_file = Path(
    "summary/functional_annotation/"
    "LCORL_missense_by_segment_haplotype.tsv"
)

# Read marker positions defining each segment
segment_positions = defaultdict(list)

with segment_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        segment_positions[row["Segment"]].append(int(row["Position"]))

for segment in segment_positions:
    segment_positions[segment] = sorted(set(segment_positions[segment]))

# Read six missense variants and their best-linked segments
missense_variants = []

with evidence_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        missense_variants.append({
            "position": int(row["Genomic_position"]),
            "protein_change": row["Protein_change"],
            "ref": row["Ref"],
            "alt": row["Alt"],
            "segment": row["Best_segment"],
        })

needed_positions = {
    variant["position"] for variant in missense_variants
}

for variant in missense_variants:
    needed_positions.update(segment_positions[variant["segment"]])

# Read phased alleles from VCF
sample_names = []
phased_alleles = {}

with gzip.open(vcf_file, "rt") as handle:
    for line in handle:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            fields = line.rstrip().split("\t")
            sample_names = fields[9:]
            continue

        fields = line.rstrip().split("\t")
        position = int(fields[1])

        if position not in needed_positions:
            continue

        ref = fields[3]
        alts = fields[4].split(",")
        allele_lookup = [ref] + alts

        format_fields = fields[8].split(":")
        gt_index = format_fields.index("GT")

        sample_pairs = []

        for sample_field in fields[9:]:
            values = sample_field.split(":")
            genotype = values[gt_index]

            # Require phased genotypes
            if "|" not in genotype:
                sample_pairs.append((None, None))
                continue

            codes = genotype.split("|")
            pair = []

            for code in codes:
                if code == ".":
                    pair.append(None)
                else:
                    index = int(code)
                    pair.append(
                        allele_lookup[index]
                        if index < len(allele_lookup)
                        else None
                    )

            sample_pairs.append(tuple(pair))

        phased_alleles[position] = sample_pairs

# Confirm every required position was found
missing_positions = sorted(
    needed_positions.difference(phased_alleles)
)

if missing_positions:
    raise RuntimeError(
        "Positions absent from phased VCF: " +
        ",".join(map(str, missing_positions))
    )

results = []

for variant in missense_variants:
    segment = variant["segment"]
    missense_position = variant["position"]
    marker_positions = segment_positions[segment]

    counts = defaultdict(
        lambda: {
            "total": 0,
            "ref": 0,
            "alt": 0,
            "other": 0,
        }
    )

    for sample_index in range(len(sample_names)):
        for chromosome_copy in (0, 1):
            marker_alleles = [
                phased_alleles[position][sample_index][chromosome_copy]
                for position in marker_positions
            ]

            missense_allele = phased_alleles[
                missense_position
            ][sample_index][chromosome_copy]

            if (
                any(allele is None for allele in marker_alleles)
                or missense_allele is None
            ):
                continue

            haplotype = "-".join(marker_alleles)
            counts[haplotype]["total"] += 1

            if missense_allele == variant["ref"]:
                counts[haplotype]["ref"] += 1
            elif missense_allele == variant["alt"]:
                counts[haplotype]["alt"] += 1
            else:
                counts[haplotype]["other"] += 1

    for haplotype, count_data in counts.items():
        total = count_data["total"]
        alt_frequency = (
            count_data["alt"] / total if total else 0
        )

        results.append({
            "Segment": segment,
            "Missense_position": missense_position,
            "Protein_change": variant["protein_change"],
            "Segment_haplotype": haplotype,
            "Haplotype_count": total,
            "Ref": variant["ref"],
            "Alt": variant["alt"],
            "Ref_count": count_data["ref"],
            "Alt_count": count_data["alt"],
            "Other_count": count_data["other"],
            "ALT_frequency_within_haplotype":
                f"{alt_frequency:.6f}",
        })

def segment_number(segment):
    return int(segment.lstrip("S"))

results.sort(
    key=lambda row: (
        segment_number(row["Segment"]),
        int(row["Missense_position"]),
        -int(row["Haplotype_count"]),
        row["Segment_haplotype"],
    )
)

columns = [
    "Segment",
    "Missense_position",
    "Protein_change",
    "Segment_haplotype",
    "Haplotype_count",
    "Ref",
    "Alt",
    "Ref_count",
    "Alt_count",
    "Other_count",
    "ALT_frequency_within_haplotype",
]

with output_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(results)

print(f"Created: {output_file}")
print(f"Rows: {len(results)}")

for variant in missense_variants:
    relevant = [
        row for row in results
        if int(row["Missense_position"]) == variant["position"]
    ]

    total = sum(int(row["Haplotype_count"]) for row in relevant)
    alt_total = sum(int(row["Alt_count"]) for row in relevant)

    print(
        variant["segment"],
        variant["protein_change"],
        f"chromosome copies={total}",
        f"ALT copies={alt_total}",
    )
