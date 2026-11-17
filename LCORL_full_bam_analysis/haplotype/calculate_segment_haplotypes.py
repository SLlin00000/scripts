import gzip
import os
from collections import Counter

VCF = "haplotype/LCORL_gene_phased_MAF05.vcf.gz"
SEGMENTS = "haplotype/LCORL_adjacent_highLD_segments.tsv"

OUTDIR = "haplotype/segment_haplotypes"
os.makedirs(OUTDIR, exist_ok=True)

# --------------------------------------------------
# Read segment definitions
# --------------------------------------------------

segments = {}

with open(SEGMENTS) as f:
    header = f.readline()

    for line in f:
        if not line.strip():
            continue

        p = line.rstrip().split("\t")

        seg_id = p[0]
        start = int(p[1])
        end = int(p[2])
        n_snps = int(p[3])
        span = int(p[4])
        positions = [int(x) for x in p[5].split(",")]

        segments[seg_id] = {
            "start": start,
            "end": end,
            "n_snps": n_snps,
            "span": span,
            "positions": positions
        }

wanted_positions = set()

for seg in segments.values():
    wanted_positions.update(seg["positions"])

# --------------------------------------------------
# Read phased VCF
# Store REF/ALT and phased GT for required positions
# --------------------------------------------------

samples = []
variants = {}

with gzip.open(VCF, "rt") as f:

    for line in f:

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            fields = line.rstrip().split("\t")
            samples = fields[9:]
            continue

        fields = line.rstrip().split("\t")

        pos = int(fields[1])

        if pos not in wanted_positions:
            continue

        ref = fields[3]
        alt = fields[4]

        fmt = fields[8].split(":")
        gt_idx = fmt.index("GT")

        gts = []

        for sample_field in fields[9:]:

            values = sample_field.split(":")

            if gt_idx >= len(values):
                gt = ".|."
            else:
                gt = values[gt_idx]

            gts.append(gt)

        variants[pos] = {
            "ref": ref,
            "alt": alt,
            "gts": gts
        }

# --------------------------------------------------
# Check positions were found
# --------------------------------------------------

missing_positions = sorted(wanted_positions - set(variants))

if missing_positions:
    print("WARNING: positions missing from VCF:")
    print(missing_positions)

print("Samples:", len(samples))
print("Segments:", len(segments))
print("Required SNPs:", len(wanted_positions))
print("VCF SNPs found:", len(variants))

# --------------------------------------------------
# Helper: allele index -> nucleotide
# --------------------------------------------------

def allele_base(pos, allele_code):

    if allele_code == ".":
        return "N"

    var = variants[pos]

    alleles = [var["ref"]] + var["alt"].split(",")

    try:
        idx = int(allele_code)
        return alleles[idx]
    except:
        return "N"

# --------------------------------------------------
# Per-segment haplotype counting
# --------------------------------------------------

summary_rows = []

for seg_id, seg in segments.items():

    positions = seg["positions"]

    hap_counts = Counter()

    sample_diplotypes = []

    valid_chromosomes = 0

    for sample_idx, sample in enumerate(samples):

        hap1 = []
        hap2 = []

        valid1 = True
        valid2 = True

        for pos in positions:

            if pos not in variants:
                valid1 = False
                valid2 = False
                continue

            gt = variants[pos]["gts"][sample_idx]

            if "|" not in gt:
                valid1 = False
                valid2 = False
                continue

            a1, a2 = gt.split("|")

            b1 = allele_base(pos, a1)
            b2 = allele_base(pos, a2)

            if b1 == "N":
                valid1 = False
            if b2 == "N":
                valid2 = False

            hap1.append(b1)
            hap2.append(b2)

        h1 = "-".join(hap1)
        h2 = "-".join(hap2)

        if valid1:
            hap_counts[h1] += 1
            valid_chromosomes += 1

        if valid2:
            hap_counts[h2] += 1
            valid_chromosomes += 1

        sample_diplotypes.append((sample, h1, h2))

    # ----------------------------------------------
    # Frequency table
    # ----------------------------------------------

    freq_out = os.path.join(
        OUTDIR,
        f"{seg_id}_haplotype_frequencies.tsv"
    )

    ordered = hap_counts.most_common()

    with open(freq_out, "w") as out:

        out.write(
            "SEGMENT\tRANK\tHAPLOTYPE\tCOUNT\tFREQUENCY\n"
        )

        for rank, (hap, count) in enumerate(ordered, 1):

            freq = (
                count / valid_chromosomes
                if valid_chromosomes > 0
                else 0
            )

            out.write(
                f"{seg_id}\t{rank}\t{hap}\t"
                f"{count}\t{freq:.6f}\n"
            )

    # ----------------------------------------------
    # Sample diplotypes
    # ----------------------------------------------

    dip_out = os.path.join(
        OUTDIR,
        f"{seg_id}_sample_diplotypes.tsv"
    )

    with open(dip_out, "w") as out:

        out.write(
            "SAMPLE\tHAPLOTYPE_1\tHAPLOTYPE_2\n"
        )

        for sample, h1, h2 in sample_diplotypes:
            out.write(
                f"{sample}\t{h1}\t{h2}\n"
            )

    # ----------------------------------------------
    # Segment summary
    # ----------------------------------------------

    n_haplotypes = len(ordered)

    if ordered:
        top_hap = ordered[0][0]
        top_count = ordered[0][1]
        top_freq = top_count / valid_chromosomes
    else:
        top_hap = "NA"
        top_count = 0
        top_freq = 0

    summary_rows.append(
        (
            seg_id,
            seg["start"],
            seg["end"],
            seg["n_snps"],
            seg["span"],
            valid_chromosomes,
            n_haplotypes,
            top_hap,
            top_count,
            top_freq
        )
    )

# --------------------------------------------------
# Global summary
# --------------------------------------------------

summary_file = os.path.join(
    OUTDIR,
    "segment_haplotype_summary.tsv"
)

with open(summary_file, "w") as out:

    out.write(
        "SEGMENT\tSTART\tEND\tN_SNPS\tSPAN_BP\t"
        "VALID_CHROMOSOMES\tN_HAPLOTYPES\t"
        "TOP_HAPLOTYPE\tTOP_COUNT\tTOP_FREQUENCY\n"
    )

    for row in summary_rows:

        out.write(
            "\t".join(map(str, row[:-1]))
            + f"\t{row[-1]:.6f}\n"
        )

print("DONE")
print("Summary:", summary_file)

