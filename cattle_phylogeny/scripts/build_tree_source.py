local = set(x.strip() for x in open("metadata/local174_sample_ids.txt") if x.strip())
run8  = set(x.strip() for x in open("metadata/run8_1089_sample_ids.txt") if x.strip())
tree  = [x.strip() for x in open("metadata/tree1217_sample_ids.txt") if x.strip()]

counts = {"LOCAL174": 0, "RUN8": 0, "UNKNOWN": 0}

with open("metadata/tree1217_source.tsv", "w") as out:
    out.write("sample_id\tsource\n")

    for s in tree:
        if s in local:
            source = "LOCAL174"
        elif s in run8:
            source = "RUN8"
        else:
            source = "UNKNOWN"

        counts[source] += 1
        out.write(f"{s}\t{source}\n")

print("Tree samples:", len(tree))
for k, v in counts.items():
    print(k + ":", v)
