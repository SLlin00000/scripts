from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

root = Path(__file__).resolve().parents[1]
records = list(SeqIO.parse(
    root / "alignments/LCORL_gene_12taxa.mafft.fasta", "fasta"
))
seqs = [str(r.seq).upper() for r in records]
assert len(records) == 12
assert len({r.id for r in records}) == 12
assert len({len(s) for s in seqs}) == 1

valid = set("ACGT")
coverage = [
    sum(base in valid for base in column)
    for column in zip(*seqs)
]
out = root / "trimmed"
out.mkdir(exist_ok=True)

for minimum in [10, 12]:
    positions = [
        i for i, count in enumerate(coverage)
        if count >= minimum
    ]
    assert positions, "No columns retained"

    filtered = []
    for record, seq in zip(records, seqs):
        new_seq = "".join(seq[i] for i in positions)
        filtered.append(SeqRecord(
            Seq(new_seq), id=record.id, description=""
        ))

    name = f"LCORL_gene_12taxa.occ{minimum}of12"
    SeqIO.write(filtered, out / f"{name}.fasta", "fasta")

    mapping = root / "metadata" / f"{name}.column_map.tsv"
    with mapping.open("w") as handle:
        handle.write(
            "filtered_column_1based\toriginal_column_1based\t"
            "ACGT_coverage\n"
        )
        for new_pos, old_pos in enumerate(positions, 1):
            handle.write(
                f"{new_pos}\t{old_pos+1}\t{coverage[old_pos]}\n"
            )

    variable = informative = 0
    for column in zip(*(str(r.seq) for r in filtered)):
        counts = Counter(b for b in column if b in valid)
        variable += len(counts) > 1
        informative += sum(n >= 2 for n in counts.values()) >= 2

    print(f"\n--- At least {minimum}/12 taxa ---")
    print("Retained columns:", len(positions))
    print("Variable columns:", variable)
    print("Parsimony-informative columns:", informative)

    for record in filtered:
        seq = str(record.seq)
        bases = sum(b in valid for b in seq)
        print(
            f"{record.id}: ACGT={bases:,} | "
            f"gap={100*seq.count('-')/len(seq):.2f}%"
        )

    print("FASTA:", out / f"{name}.fasta")
    print("Column map:", mapping)
