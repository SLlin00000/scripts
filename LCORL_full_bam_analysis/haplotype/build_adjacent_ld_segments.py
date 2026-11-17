VCF_POS_FILE = "haplotype/common_snp_positions.txt"
LD_FILE = "haplotype/LCORL_gene_MAF05_LD.hap.ld"
OUTFILE = "haplotype/LCORL_adjacent_highLD_segments.tsv"

# --------------------------------------------------
# Read ordered SNP positions
# --------------------------------------------------

positions = []

with open(VCF_POS_FILE) as f:
    for line in f:
        if line.strip():
            positions.append(int(line.strip()))

positions = sorted(set(positions))

# --------------------------------------------------
# Read pairwise LD
# key = (smaller_position, larger_position)
# --------------------------------------------------

ld = {}

with open(LD_FILE) as f:
    header = f.readline()

    for line in f:
        if not line.strip():
            continue

        p = line.split()

        pos1 = int(p[1])
        pos2 = int(p[2])

        try:
            r2 = float(p[4])
        except ValueError:
            continue

        if r2 != r2:   # NaN
            continue

        key = tuple(sorted((pos1, pos2)))
        ld[key] = r2

# --------------------------------------------------
# Build strict runs:
# consecutive genomic SNPs must have r2 >= 0.8
# --------------------------------------------------

segments = []
current = []

for i in range(len(positions)):

    if not current:
        current = [positions[i]]

    if i == len(positions) - 1:
        break

    p1 = positions[i]
    p2 = positions[i + 1]

    r2 = ld.get((p1, p2))

    if r2 is not None and r2 >= 0.8:
        current.append(p2)

    else:
        if len(current) >= 3:
            segments.append(current)

        current = [p2]

if len(current) >= 3:
    segments.append(current)

# --------------------------------------------------
# Output
# --------------------------------------------------

with open(OUTFILE, "w") as out:

    out.write(
        "SEGMENT\tSTART\tEND\tN_SNPS\tSPAN_BP\tSNP_POSITIONS\n"
    )

    for i, seg in enumerate(segments, 1):

        out.write(
            f"S{i}\t"
            f"{seg[0]}\t"
            f"{seg[-1]}\t"
            f"{len(seg)}\t"
            f"{seg[-1]-seg[0]}\t"
            + ",".join(map(str, seg))
            + "\n"
        )

print("Total common SNPs:", len(positions))
print("Strict adjacent high-LD segments:", len(segments))

for i, seg in enumerate(segments, 1):
    print(
        f"S{i}: "
        f"{seg[0]}-{seg[-1]} "
        f"N={len(seg)} "
        f"span={seg[-1]-seg[0]} bp"
    )

