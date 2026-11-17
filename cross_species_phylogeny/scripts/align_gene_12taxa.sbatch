#!/bin/bash
#SBATCH --job-name=LCORL_gene12
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/gene12_%j.out
#SBATCH --error=logs/gene12_%j.err

set -euo pipefail

cd /QRISdata/Q9774/lcorl_project/cross_species_phylogeny

module load mafft/7.526-gcc-13.3.0-with-extensions

mafft --version > logs/gene12_mafft_version.txt 2>&1

mafft \
    --nuc \
    --thread "$SLURM_CPUS_PER_TASK" \
    --memsave \
    --retree 2 \
    --maxiterate 0 \
    selected/LCORL_gene_12taxa.fasta \
    > alignments/LCORL_gene_12taxa.mafft.fasta.tmp \
    2> logs/gene12_mafft.log

python3 - <<'PY'
from Bio import SeqIO

source = list(SeqIO.parse("selected/LCORL_gene_12taxa.fasta", "fasta"))
aligned = list(SeqIO.parse(
    "alignments/LCORL_gene_12taxa.mafft.fasta.tmp", "fasta"
))

expected = {r.id: str(r.seq).upper() for r in source}
observed = {r.id: str(r.seq).upper().replace("-", "") for r in aligned}

assert len(aligned) == 12, "Incorrect sequence count"
assert len(observed) == 12, "Duplicate sequence IDs"
assert expected == observed, "Sequences changed during alignment"
assert len({len(r.seq) for r in aligned}) == 1, "Unequal alignment lengths"

print("PASS: 12 aligned sequences; original bases preserved.")
print("Alignment length:", len(aligned[0].seq))
PY

mv alignments/LCORL_gene_12taxa.mafft.fasta.tmp \
   alignments/LCORL_gene_12taxa.mafft.fasta
