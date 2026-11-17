import os

INFILE = "haplotype/LCORL_gene_MAF05_highLD_r2_08.tsv"
OUTDIR = "haplotype/candidate_blocks"
os.makedirs(OUTDIR, exist_ok=True)

class UnionFind:
    def __init__(self):
        self.parent = {}
    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

pairs = []
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
            pairs.append((pos1, pos2, r2))

for max_dist in [20000, 50000, 100000]:
    uf = UnionFind()
    used_pairs = 0

    for p1, p2, r2 in pairs:
        if abs(p2 - p1) <= max_dist:
            uf.union(p1, p2)
            used_pairs += 1

    groups = {}
    for snp in uf.parent:
        root = uf.find(snp)
        groups.setdefault(root, set()).add(snp)

    blocks = []
    for snps in groups.values():
        if len(snps) >= 3:
            s = sorted(snps)
            blocks.append((s[0], s[-1], len(s), s[-1]-s[0], s))

    blocks.sort()

    summary = os.path.join(
        OUTDIR,
        f"candidate_blocks_r2_08_maxdist_{max_dist//1000}kb.tsv"
    )

    with open(summary, "w") as out:
        out.write("BLOCK\tSTART\tEND\tN_SNPS\tSPAN_BP\tSNP_POSITIONS\n")
        for i, (start, end, n, span, snps) in enumerate(blocks, 1):
            out.write(
                f"B{i}\t{start}\t{end}\t{n}\t{span}\t"
                + ",".join(map(str, snps)) + "\n"
            )

    print(
        f"max_dist={max_dist//1000}kb",
        f"pairs_used={used_pairs}",
        f"candidate_blocks={len(blocks)}"
    )

print("DONE")
