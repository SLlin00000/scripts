import os

INFILE = "haplotype/LCORL_gene_MAF05_highLD_r2_08.tsv"
OUTFILE = "haplotype/LCORL_candidate_segments_20kb.tsv"

positions = set()

with open(INFILE) as f:
    header = f.readline()
    for line in f:
        if not line.strip():
            continue
        p = line.split()
        pos1 = int(p[1])
        pos2 = int(p[2])
        r2 = float(p[4])

        if r2 >= 0.8:
            positions.add(pos1)
            positions.add(pos2)

positions = sorted(positions)

segments = []
current = []

for pos in positions:
    if not current:
        current = [pos]
        continue

    if pos - current[-1] <= 20000:
        current.append(pos)
    else:
        if len(current) >= 3:
            segments.append(current)
        current = [pos]

if current and len(current) >= 3:
    segments.append(current)

with open(OUTFILE, "w") as out:
    out.write("SEGMENT\tSTART\tEND\tN_SNPS\tSPAN_BP\tSNP_POSITIONS\n")

    for i, seg in enumerate(segments, 1):
        out.write(
            f"S{i}\t{seg[0]}\t{seg[-1]}\t{len(seg)}\t"
            f"{seg[-1]-seg[0]}\t"
            + ",".join(map(str, seg))
            + "\n"
        )

print("Segments:", len(segments))
print("Output:", OUTFILE)

