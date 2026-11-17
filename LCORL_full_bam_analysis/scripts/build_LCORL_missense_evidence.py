#!/usr/bin/env python3

import csv
import gzip
import re
from pathlib import Path

vep_file = Path("vep/LCORL_representative_missense_VEP107.tsv")
ld_file = Path("haplotype/LCORL_missense_segment_maxLD.tsv")
vcf_file = Path("haplotype/LCORL_common_CDS_snps.vcf.gz")

out_file = Path(
    "summary/functional_annotation/"
    "LCORL_missense_integrated_evidence.tsv"
)

out_file.parent.mkdir(parents=True, exist_ok=True)

# 1. Read VEP results
variants = {}

with vep_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        position = int(row["Location"].split(":")[-1].split("-")[0])
        protein_position = row["Protein_position"]
        ref_aa, alt_aa = row["Amino_acids"].split("/")

        sift_match = re.match(r"(.+)\(([0-9.]+)\)", row["SIFT"])
        if sift_match:
            sift_prediction = sift_match.group(1)
            sift_score = float(sift_match.group(2))
        else:
            sift_prediction = row["SIFT"]
            sift_score = None

        variants[position] = {
            "Genomic_position": position,
            "Protein_change": f"{ref_aa}{protein_position}{alt_aa}",
            "SIFT_prediction": sift_prediction,
            "SIFT_score": sift_score,
            "DUF4553_status": (
                "Within_DUF4553"
                if "PF15090" in row["Domains"]
                else "Outside_DUF4553"
            ),
            "VEP_domains": row["Domains"],
        }

# 2. Add haplotype/LD evidence
with ld_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        position = int(row["MISSENSE_POS"])

        if position in variants:
            variants[position].update({
                "Best_segment": row["BEST_SEGMENT"],
                "Best_segment_SNP": row["BEST_SEGMENT_SNP"],
                "Max_r2": row["MAX_R2"],
                "D": row["D"],
                "Dprime": row["Dprime"],
                "LD_class": row["LD_CLASS"],
            })

# 3. Calculate AC, AN, AF and MAF from VCF genotypes
with gzip.open(vcf_file, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")
        position = int(fields[1])

        if position not in variants:
            continue

        ref = fields[3]
        alt = fields[4]

        format_fields = fields[8].split(":")
        gt_index = format_fields.index("GT")

        alt_count = 0
        allele_number = 0
        called_samples = 0

        for sample in fields[9:]:
            sample_fields = sample.split(":")

            if gt_index >= len(sample_fields):
                continue

            genotype = sample_fields[gt_index]

            if genotype in {".", "./.", ".|."}:
                continue

            alleles = genotype.replace("|", "/").split("/")
            valid_alleles = [a for a in alleles if a != "."]

            if not valid_alleles:
                continue

            called_samples += 1

            for allele in valid_alleles:
                allele_number += 1
                if int(allele) > 0:
                    alt_count += 1

        allele_frequency = (
            alt_count / allele_number if allele_number else None
        )
        maf = (
            min(allele_frequency, 1 - allele_frequency)
            if allele_frequency is not None
            else None
        )

        variants[position].update({
            "Ref": ref,
            "Alt": alt,
            "Called_samples": called_samples,
            "AC": alt_count,
            "AN": allele_number,
            "AF": allele_frequency,
            "MAF": maf,
        })

# 4. Write integrated table
columns = [
    "Genomic_position",
    "Ref",
    "Alt",
    "Protein_change",
    "SIFT_prediction",
    "SIFT_score",
    "DUF4553_status",
    "Best_segment",
    "Best_segment_SNP",
    "Max_r2",
    "D",
    "Dprime",
    "LD_class",
    "Called_samples",
    "AC",
    "AN",
    "AF",
    "MAF",
    "VEP_domains",
]

with out_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter="\t",
        extrasaction="ignore",
    )
    writer.writeheader()

    for position in sorted(variants):
        row = variants[position]

        for field in ("AF", "MAF"):
            if isinstance(row.get(field), float):
                row[field] = f"{row[field]:.6f}"

        writer.writerow(row)

print(f"Created: {out_file}")
print(f"Variants: {len(variants)}")

