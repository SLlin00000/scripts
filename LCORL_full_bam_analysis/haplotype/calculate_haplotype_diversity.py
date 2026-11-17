import glob
import os

indir = "haplotype/segment_haplotypes"
outfile = f"{indir}/segment_haplotype_diversity.tsv"

results = []

for path in sorted(glob.glob(f"{indir}/S*_haplotype_frequencies.tsv")):

    segment = os.path.basename(path).split("_")[0]

    freqs = []

    with open(path) as f:
        header = f.readline()

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")
            freq = float(fields[4])
            freqs.append(freq)

    if not freqs:
        continue

    diversity = 1 - sum(p*p for p in freqs)

    effective_haplotypes = (
        1 / sum(p*p for p in freqs)
    )

    top_freq = max(freqs)

    common_haps = sum(
        1 for p in freqs if p >= 0.05
    )

    results.append(
        (
            segment,
            len(freqs),
            common_haps,
            top_freq,
            diversity,
            effective_haplotypes
        )
    )

def seg_num(x):
    return int(x[0][1:])

results.sort(key=seg_num)

with open(outfile, "w") as out:

    out.write(
        "SEGMENT\tN_HAPLOTYPES\tN_HAPLOTYPES_FREQ_GE_0.05\t"
        "TOP_FREQUENCY\tHAPLOTYPE_DIVERSITY\t"
        "EFFECTIVE_N_HAPLOTYPES\n"
    )

    for row in results:

        out.write(
            f"{row[0]}\t"
            f"{row[1]}\t"
            f"{row[2]}\t"
            f"{row[3]:.6f}\t"
            f"{row[4]:.6f}\t"
            f"{row[5]:.3f}\n"
        )

print("DONE")
print(outfile)

