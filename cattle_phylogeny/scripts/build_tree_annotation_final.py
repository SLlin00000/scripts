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

        ancestry = row["expected_ancestry"].strip() or "metadata_unavailable"
        production = row["production_type"].strip() or "metadata_unavailable"
        breed = row["normalised_breed"].strip() or "metadata_unavailable"

        run8_meta[sid] = {
            "ancestry": ancestry,
            "production": production,
            "breed": breed
        }

with open("metadata/tree1217_annotation_final.tsv", "w") as out:
    out.write(
        "sample_id\tsource\tancestry_group\tproduction_type\tbreed\t"
        "annotation_status\tlabel_flag\n"
    )

    for sid in tree_ids:
        if sid in local:
            source = "LOCAL174"
            ancestry = "metadata_unavailable"
            production = "metadata_unavailable"
            breed = "metadata_unavailable"
            annotation_status = "source_only"

        elif sid in run8_meta:
            source = "RUN8"
            ancestry = run8_meta[sid]["ancestry"]
            production = run8_meta[sid]["production"]
            breed = run8_meta[sid]["breed"]
            annotation_status = "annotated"

        else:
            source = "UNKNOWN"
            ancestry = "metadata_unavailable"
            production = "metadata_unavailable"
            breed = "metadata_unavailable"
            annotation_status = "unmatched"

        out.write(
            f"{sid}\t{source}\t{ancestry}\t{production}\t{breed}\t"
            f"{annotation_status}\tFALSE\n"
        )

print("Written: metadata/tree1217_annotation_final.tsv")
