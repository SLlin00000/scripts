segments_file = "haplotype/LCORL_adjacent_highLD_segments.tsv"
ld_file = "haplotype/LCORL_missense_vs_segments_LD.tsv"
out_file = "haplotype/LCORL_missense_segment_maxLD.tsv"

# segment SNP -> segment ID
snp_to_segment = {}

with open(segments_file) as f:
    header = f.readline()

    for line in f:
        if not line.strip():
            continue

        parts = line.rstrip().split("\t")
        segment = parts[0]

        for pos in parts[5].split(","):
            snp_to_segment[int(pos)] = segment


best = {}

with open(ld_file) as f:
    header = f.readline()

    for line in f:
        if not line.strip():
            continue

        missense_pos, segment_pos, r2, d, dprime = line.rstrip().split("\t")

        missense_pos = int(missense_pos)
        segment_pos = int(segment_pos)
        r2 = float(r2)

        segment = snp_to_segment.get(segment_pos, "NA")

        if missense_pos not in best or r2 > best[missense_pos]["r2"]:
            best[missense_pos] = {
                "segment": segment,
                "segment_pos": segment_pos,
                "r2": r2,
                "d": d,
                "dprime": dprime
            }


with open(out_file, "w") as out:

    out.write(
        "MISSENSE_POS\tBEST_SEGMENT\tBEST_SEGMENT_SNP\t"
        "MAX_R2\tD\tDprime\tLD_CLASS\n"
    )

    for pos in sorted(best):

        x = best[pos]

        if x["r2"] >= 0.8:
            ld_class = "strong"
        elif x["r2"] >= 0.5:
            ld_class = "moderate"
        else:
            ld_class = "weak"

        out.write(
            f"{pos}\t{x['segment']}\t{x['segment_pos']}\t"
            f"{x['r2']:.6f}\t{x['d']}\t{x['dprime']}\t"
            f"{ld_class}\n"
        )

print("DONE")
print(out_file)
