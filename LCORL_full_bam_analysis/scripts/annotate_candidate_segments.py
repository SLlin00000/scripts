import gzip
from collections import defaultdict

SEGMENTS = "haplotype/LCORL_adjacent_highLD_segments.tsv"
GFF = "annotation/Bos_taurus.ARS-UCD1.2.105.chr.gff3.gz"
OUT = "haplotype/segment_SNP_annotation.tsv"

LCORL_START = 37377806
LCORL_END = 37556970

# --------------------------------------------------
# Read SNP positions from S1-S12
# --------------------------------------------------

snps = []

with open(SEGMENTS) as f:
    header = f.readline()

    for line in f:
        if not line.strip():
            continue

        parts = line.rstrip().split("\t")

        segment = parts[0]
        positions = parts[5].split(",")

        for p in positions:
            snps.append((segment, int(p)))

# --------------------------------------------------
# Read Ensembl features overlapping LCORL
# --------------------------------------------------

features = defaultdict(list)

with gzip.open(GFF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.rstrip().split("\t")

        if len(parts) < 9:
            continue

        chrom, source, feature, start, end, score, strand, phase, attrs = parts

        if chrom != "6":
            continue

        start = int(start)
        end = int(end)

        if end < LCORL_START or start > LCORL_END:
            continue

        if feature not in [
            "exon",
            "CDS",
            "five_prime_UTR",
            "three_prime_UTR"
        ]:
            continue

        features[feature].append((start, end, strand, attrs))

# --------------------------------------------------
# Annotate SNPs
# --------------------------------------------------

def overlapping(feature_type, pos):
    hits = []

    for start, end, strand, attrs in features[feature_type]:
        if start <= pos <= end:
            hits.append((start, end, strand, attrs))

    return hits


with open(OUT, "w") as out:

    out.write(
        "Segment\tPosition\tFeature\tCDS\t"
        "Exon_hits\tCDS_hits\tUTR_hits\n"
    )

    for segment, pos in snps:

        cds_hits = overlapping("CDS", pos)
        exon_hits = overlapping("exon", pos)

        utr5_hits = overlapping("five_prime_UTR", pos)
        utr3_hits = overlapping("three_prime_UTR", pos)

        utr_hits = utr5_hits + utr3_hits

        if cds_hits:
            feature = "CDS"
            cds = "YES"

        elif utr_hits:
            feature = "UTR"
            cds = "NO"

        elif exon_hits:
            feature = "exon_nonCDS"
            cds = "NO"

        else:
            feature = "intron"
            cds = "NO"

        out.write(
            f"{segment}\t{pos}\t{feature}\t{cds}\t"
            f"{len(exon_hits)}\t{len(cds_hits)}\t{len(utr_hits)}\n"
        )

print("DONE")
print(OUT)

