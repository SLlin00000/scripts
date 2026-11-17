#!/usr/bin/env python3

import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

vcf_file = Path("haplotype/LCORL_gene_phased_MAF01.vcf.gz")
evidence_file = Path(
    "summary/functional_annotation/"
    "LCORL_representative_missense_MAF01_evidence.tsv"
)
segment_file = Path("haplotype/segment_SNP_annotation.tsv")

protein_output = Path(
    "summary/functional_annotation/"
    "LCORL_joint_missense_protein_haplotypes_MAF01.tsv"
)
mapping_output = Path(
    "summary/functional_annotation/"
    "LCORL_joint_missense_segment_mapping_MAF01.tsv"
)

segments_to_use = ["S1", "S5", "S11"]

# Read missense variants
variants = []

with evidence_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        match = re.fullmatch(
            r"([A-Z*])(\d+)([A-Z*])",
            row["Protein_change"],
        )

        if not match:
            raise ValueError(
                f"Cannot parse protein change: {row['Protein_change']}"
            )

        variants.append({
            "position": int(row["Genomic_position"]),
            "ref_nt": row["Ref"],
            "alt_nt": row["Alt"],
            "ref_aa": match.group(1),
            "protein_position": int(match.group(2)),
            "alt_aa": match.group(3),
            "change": row["Protein_change"],
        })

variants.sort(key=lambda item: item["protein_position"])

# Read marker positions for S1, S5 and S11
segment_positions = defaultdict(list)

with segment_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        segment = row["Segment"]

        if segment in segments_to_use:
            segment_positions[segment].append(int(row["Position"]))

for segment in segment_positions:
    segment_positions[segment] = sorted(
        set(segment_positions[segment])
    )

needed_positions = {item["position"] for item in variants}

for segment in segments_to_use:
    needed_positions.update(segment_positions[segment])

# Read phased VCF
samples = []
alleles_by_position = {}

with gzip.open(vcf_file, "rt") as handle:
    for line in handle:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            fields = line.rstrip().split("\t")
            samples = fields[9:]
            continue

        fields = line.rstrip().split("\t")
        position = int(fields[1])

        if position not in needed_positions:
            continue

        allele_lookup = [fields[3]] + fields[4].split(",")
        format_fields = fields[8].split(":")
        gt_index = format_fields.index("GT")

        sample_pairs = []

        for sample_field in fields[9:]:
            genotype = sample_field.split(":")[gt_index]

            if "|" not in genotype:
                sample_pairs.append((None, None))
                continue

            pair = []

            for code in genotype.split("|"):
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

        alleles_by_position[position] = sample_pairs

missing = sorted(needed_positions.difference(alleles_by_position))

if missing:
    raise RuntimeError(
        "Positions missing from VCF: " +
        ",".join(map(str, missing))
    )

protein_counts = Counter()
mapping_counts = Counter()
usable_chromosomes = 0

for sample_index in range(len(samples)):
    for chromosome_copy in (0, 1):
        amino_acid_states = []
        alt_changes = []
        valid = True

        for variant in variants:
            nucleotide = alleles_by_position[
                variant["position"]
            ][sample_index][chromosome_copy]

            if nucleotide == variant["ref_nt"]:
                amino_acid = variant["ref_aa"]
            elif nucleotide == variant["alt_nt"]:
                amino_acid = variant["alt_aa"]
                alt_changes.append(variant["change"])
            else:
                valid = False
                break

            amino_acid_states.append(
                f"{variant['protein_position']}{amino_acid}"
            )

        if not valid:
            continue

        segment_haplotypes = {}

        for segment in segments_to_use:
            marker_alleles = [
                alleles_by_position[position][
                    sample_index
                ][chromosome_copy]
                for position in segment_positions[segment]
            ]

            if any(allele is None for allele in marker_alleles):
                valid = False
                break

            segment_haplotypes[segment] = "-".join(marker_alleles)

        if not valid:
            continue

        protein_haplotype = "-".join(amino_acid_states)
        change_pattern = (
            ",".join(alt_changes)
            if alt_changes
            else f"Reference_at_all_{len(variants)}"
        )

        protein_counts[
            (protein_haplotype, change_pattern)
        ] += 1

        mapping_counts[(
            protein_haplotype,
            change_pattern,
            segment_haplotypes["S1"],
            segment_haplotypes["S5"],
            segment_haplotypes["S11"],
        )] += 1

        usable_chromosomes += 1

# Write protein haplotype frequencies
with protein_output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow([
        "Rank",
        "Amino_acid_haplotype",
        "ALT_missense_changes",
        "Count",
        "Frequency",
    ])

    ranked = sorted(
        protein_counts.items(),
        key=lambda item: (-item[1], item[0]),
    )

    for rank, ((protein_haplotype, changes), count) in enumerate(
        ranked,
        start=1,
    ):
        writer.writerow([
            rank,
            protein_haplotype,
            changes,
            count,
            f"{count / usable_chromosomes:.6f}",
        ])

# Write mapping to segment haplotypes
with mapping_output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow([
        "Amino_acid_haplotype",
        "ALT_missense_changes",
        "S1_haplotype",
        "S5_haplotype",
        "S11_haplotype",
        "Count",
        "Frequency",
    ])

    ranked_mapping = sorted(
        mapping_counts.items(),
        key=lambda item: (-item[1], item[0]),
    )

    for key, count in ranked_mapping:
        writer.writerow([
            *key,
            count,
            f"{count / usable_chromosomes:.6f}",
        ])

print(f"Complete chromosome copies: {usable_chromosomes}")
print(f"Distinct protein haplotypes: {len(protein_counts)}")
print(f"Created: {protein_output}")
print(f"Created: {mapping_output}")

