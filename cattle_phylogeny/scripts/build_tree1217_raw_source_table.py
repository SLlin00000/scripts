#!/usr/bin/env python3

import csv
from collections import defaultdict

TREE_IDS = "metadata/tree1217_sample_ids.txt"
RUN8_RAW = "metadata/run8_metadata_raw.tsv"
LOCAL_MAP = "metadata/q9827_bam_to_sample.tsv"

OUT = "metadata/tree1217_raw_source_table.tsv"

tree_ids = []

with open(TREE_IDS) as f:
    for line in f:
        x = line.strip()
        if x:
            tree_ids.append(x)

tree_set = set(tree_ids)

print("Tree sample IDs:", len(tree_ids))
print("Unique tree sample IDs:", len(tree_set))

run8_rows = []

with open(RUN8_RAW, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    run8_fields = reader.fieldnames

    for row in reader:
        run8_rows.append(row)

print("Run8 raw fields:")
print(run8_fields)

run8_index = {}

for row in run8_rows:
    for field in ["ENA_accession", "BioSample", "sample_name"]:
        value = row.get(field, "").strip()
        if value:
            run8_index[value] = row

local_index = defaultdict(list)

with open(LOCAL_MAP, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        sid = row["sample_id"].strip()
        if sid:
            local_index[sid].append(row)

out_fields = [
    "tree_sample_id",
    "source",
    "ENA_accession",
    "BioSample",
    "scientific_name",
    "organism_attribute",
    "isolate",
    "breed",
    "sex",
    "country",
    "sample_name",
    "dev_stage",
    "description",
    "bam_file",
    "sample_id"
]

n_run8 = 0
n_local_samples = 0
n_local_rows = 0
unmatched = []

with open(OUT, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=out_fields
    )
    writer.writeheader()

    for sid in tree_ids:

        if sid in local_index:
            n_local_samples += 1

            for local_row in local_index[sid]:
                out = {k: "" for k in out_fields}

                out["tree_sample_id"] = sid
                out["source"] = "LOCAL174"
                out["bam_file"] = local_row["bam_file"]
                out["sample_id"] = local_row["sample_id"]

                writer.writerow(out)
                n_local_rows += 1

            continue

        if sid in run8_index:
            raw = run8_index[sid]

            out = {k: "" for k in out_fields}

            out["tree_sample_id"] = sid
            out["source"] = "RUN8"

            for field in [
                "ENA_accession",
                "BioSample",
                "scientific_name",
                "organism_attribute",
                "isolate",
                "breed",
                "sex",
                "country",
                "sample_name",
                "dev_stage",
                "description"
            ]:
                out[field] = raw.get(field, "")

            writer.writerow(out)
            n_run8 += 1
            continue

        unmatched.append(sid)

print()
print("Written:", OUT)
print("RUN8 tree samples:", n_run8)
print("LOCAL tree samples:", n_local_samples)
print("LOCAL BAM rows written:", n_local_rows)
print("Unmatched tree IDs:", len(unmatched))

if unmatched:
    print("\nUnmatched IDs:")
    for x in unmatched:
        print(x)
