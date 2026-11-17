input_file = "haplotype/LCORL_missense_protein_consequences.tsv"
detail_out = "haplotype/LCORL_missense_protein_detail.tsv"
summary_out = "haplotype/LCORL_missense_protein_summary.tsv"

detail = []
by_snp = {}

with open(input_file) as f:
    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip().split("\t")

        pos = fields[0]
        ref = fields[1]
        alt = fields[2]
        bcsq = fields[3]

        by_snp.setdefault(pos, {
            "ref": ref,
            "alt": alt,
            "transcripts": set(),
            "aa_changes": set()
        })

        for consequence in bcsq.split(","):

            p = consequence.split("|")

            if len(p) < 6:
                continue

            if p[0] != "missense":
                continue

            gene = p[1]
            transcript = p[2]
            biotype = p[3]
            strand = p[4]
            aa_change = p[5]

            detail.append([
                pos,
                ref,
                alt,
                gene,
                transcript,
                biotype,
                strand,
                aa_change
            ])

            by_snp[pos]["transcripts"].add(transcript)
            by_snp[pos]["aa_changes"].add(aa_change)


with open(detail_out, "w") as out:
    out.write(
        "POS\tREF\tALT\tGENE\tTRANSCRIPT\tBIOTYPE\tSTRAND\tAA_CHANGE\n"
    )

    for row in detail:
        out.write("\t".join(row) + "\n")


with open(summary_out, "w") as out:
    out.write(
        "POS\tREF\tALT\tN_TRANSCRIPTS\tN_UNIQUE_AA_CHANGES\tAA_CHANGES\n"
    )

    for pos in sorted(by_snp, key=int):

        d = by_snp[pos]

        aa = sorted(d["aa_changes"])
        tx = sorted(d["transcripts"])

        out.write(
            f"{pos}\t{d['ref']}\t{d['alt']}\t"
            f"{len(tx)}\t{len(aa)}\t"
            f"{','.join(aa)}\n"
        )

print("DONE")
print(detail_out)
print(summary_out)

