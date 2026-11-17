import csv

tree_ids = [
    x.strip()
    for x in open("metadata/tree1217_sample_ids.txt")
    if x.strip()
]

local = set(
    x.strip()
    for x in open("metadata/local174_sample_ids.txt")
    if x.strip()
)

run8_meta = {}

with open("metadata/run8_samples_final_classification.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        sid = row["BioSample"].strip()

        ancestry = row["expected_ancestry"].strip()
        production = row["production_type"].strip()
        breed = row["normalised_breed"].strip()
        status = row["final_status"].strip()

        run8_meta[sid] = {
            "ancestry": ancestry if ancestry else "UNKNOWN",
            "production": production if production else "UNKNOWN",
            "breed": breed if breed else "UNKNOWN",
            "status": status if status else "UNKNOWN"
        }

with open("metadata/tree1217_annotation_stage1.tsv", "w") as out:
    out.write(
        "sample_id\tsource\tancestry_group\tproduction_type\t"
        "breed\tlabel_flag\n"
    )

    for sid in tree_ids:

        if sid in local:
            source = "LOCAL174"
            ancestry = "UNKNOWN"
            production = "UNKNOWN"
            breed = "UNKNOWN"

        elif sid in run8_meta:
            source = "RUN8"
            ancestry = run8_meta[sid]["ancestry"]
            production = run8_meta[sid]["production"]
            breed = run8_meta[sid]["breed"]

        else:
            source = "UNKNOWN"
            ancestry = "UNKNOWN"
            production = "UNKNOWN"
            breed = "UNKNOWN"

        out.write(
            f"{sid}\t{source}\t{ancestry}\t{production}\t"
            f"{breed}\tFALSE\n"
        )

print("Written: metadata/tree1217_annotation_stage1.tsv")
