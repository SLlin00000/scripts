#!/usr/bin/env python3

import csv

OLD = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/breed_dictionary_v1_backup.tsv"
NEW = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/current_125_breeds.tsv"
OUT = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/breed_dictionary_master.tsv"

old_map = {}

with open(OLD, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        old_map[row["raw_breed"]] = row

current = {}

with open(NEW, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        current[row["raw_breed"]] = row["count_current"]

fields = [
    "raw_breed",
    "normalised_breed",
    "expected_ancestry",
    "production_type",
    "status",
    "note",
    "count_current",
    "present_in_current"
]

with open(OUT, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()

    # current 125 first
    for breed, count in current.items():

        if breed in old_map:
            old = old_map[breed]

            writer.writerow({
                "raw_breed": breed,
                "normalised_breed": old.get("normalised_breed", breed),
                "expected_ancestry": old.get("expected_ancestry", ""),
                "production_type": old.get("production_type", ""),
                "status": old.get("status", "REVIEW"),
                "note": old.get("note", ""),
                "count_current": count,
                "present_in_current": "yes"
            })

        else:
            writer.writerow({
                "raw_breed": breed,
                "normalised_breed": breed,
                "expected_ancestry": "",
                "production_type": "",
                "status": "REVIEW",
                "note": "new label from unified raw-metadata reclassification",
                "count_current": count,
                "present_in_current": "yes"
            })

    # retain old labels no longer present
    for breed, old in old_map.items():

        if breed not in current:
            writer.writerow({
                "raw_breed": breed,
                "normalised_breed": old.get("normalised_breed", breed),
                "expected_ancestry": old.get("expected_ancestry", ""),
                "production_type": old.get("production_type", ""),
                "status": old.get("status", "REVIEW"),
                "note": (old.get("note", "") + "; not present in current reclassified set").strip("; "),
                "count_current": "0",
                "present_in_current": "no"
            })

