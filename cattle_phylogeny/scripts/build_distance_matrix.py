psam = "results/cattle1217_tree_input.psam"
pairs = "results/cattle1217_pairwise_distance.tsv"
out = "results/cattle1217_distance_matrix.tsv"

ids = []

with open(psam) as f:
    header = f.readline().strip().split()
    iid_col = header.index("#IID")

    for line in f:
        if line.strip():
            fields = line.split()
            ids.append(fields[iid_col])

index = {sample: i for i, sample in enumerate(ids)}

n = len(ids)
D = [[0.0] * n for _ in range(n)]

with open(pairs) as f:
    next(f)

    for line in f:
        iid1, iid2, dist, nsnp = line.rstrip().split("\t")

        i = index[iid1]
        j = index[iid2]
        d = float(dist)

        D[i][j] = d
        D[j][i] = d

symmetric = True

for i in range(n):
    for j in range(i + 1, n):
        if abs(D[i][j] - D[j][i]) > 1e-12:
            symmetric = False
            break

diagonal_zero = all(abs(D[i][i]) < 1e-12 for i in range(n))

with open(out, "w") as f:
    f.write("IID\t" + "\t".join(ids) + "\n")

    for i, sample in enumerate(ids):
        values = "\t".join(f"{x:.8f}" for x in D[i])
        f.write(sample + "\t" + values + "\n")

print("Samples:", n)
print("Matrix shape:", f"({n}, {n})")
print("Symmetric:", symmetric)
print("Diagonal zero:", diagonal_zero)
print("Written:", out)

