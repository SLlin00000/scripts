missense_file = "haplotype/LCORL_missense_positions.txt"
segment_file = "haplotype/LCORL_segment_positions.txt"
ld_file = "haplotype/LCORL_gene_MAF05_LD.hap.ld"
out_file = "haplotype/LCORL_missense_vs_segments_LD.tsv"

with open(missense_file) as f:
    missense = set(int(x.strip()) for x in f if x.strip())

with open(segment_file) as f:
    segments = set(int(x.strip()) for x in f if x.strip())

with open(ld_file) as f, open(out_file, "w") as out:
    header = f.readline().strip()
    out.write("MISSENSE_POS\tSEGMENT_POS\tR2\tD\tDprime\n")

    for line in f:
        parts = line.strip().split()
        if len(parts) < 7:
            continue

        pos1 = int(parts[1])
        pos2 = int(parts[2])
        r2 = parts[4]
        d = parts[5]
        dprime = parts[6]

        if pos1 in missense and pos2 in segments:
            out.write(f"{pos1}\t{pos2}\t{r2}\t{d}\t{dprime}\n")

        elif pos2 in missense and pos1 in segments:
            out.write(f"{pos2}\t{pos1}\t{r2}\t{d}\t{dprime}\n")

print("DONE")
print(out_file)

