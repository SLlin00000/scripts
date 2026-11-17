#!/usr/bin/env python3

import csv
import re

INFILE = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/run8_metadata_raw.tsv"
OUTFILE = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/run8_metadata_preclassified.tsv"

MISSING = {
    "", "na", "n/a", "none", "null",
    "missing", "unknown",
    "not provided", "not collected"
}

def clean(x):
    return (x or "").strip()

def valid(x):
    return clean(x).lower() not in MISSING

def ancestry_group(scientific_name):
    s = clean(scientific_name).lower()

    if ("bos indicus" in s and "bos taurus" in s):
        return "admixed"

    if "bos indicus" in s:
        return "indicine"

    if "bos taurus" in s:
        return "taurine"

    if "bos primigenius" in s:
        return "outgroup_or_review"

    return "unknown"


def choose_breed(row):

    # 1. Explicit breed field = strongest evidence
    if valid(row["breed"]):
        return clean(row["breed"]), "breed", "high"

    # 2. Use isolate only if explicit breed is absent
    if valid(row["isolate"]):
        return clean(row["isolate"]), "isolate", "review"

    # 3. Weak fallback
    if valid(row["sample_name"]):
        return clean(row["sample_name"]), "sample_name", "low"

    return "", "none", "unresolved"


with open(INFILE, newline="") as inp, \
     open(OUTFILE, "w", newline="") as out:

    reader = csv.DictReader(inp, delimiter="\t")

    fields = reader.fieldnames + [
        "ancestry_group",
        "breed_candidate",
        "breed_source",
        "breed_confidence",
        "label_status",
        "production_type"
    ]

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for row in reader:

        ancestry = ancestry_group(
            row["scientific_name"]
        )

        breed, source, confidence = choose_breed(row)

        if ancestry == "unknown":
            status = "needs_review"

        elif confidence == "high":
            status = "well_labelled"

        elif confidence == "review":
            status = "needs_breed_review"

        else:
            status = "poorly_labelled"

        row["ancestry_group"] = ancestry
        row["breed_candidate"] = breed
        row["breed_source"] = source
        row["breed_confidence"] = confidence
        row["label_status"] = status

        # intentionally blank:
        # beef/dairy classification will be manually verified
        row["production_type"] = ""

        writer.writerow(row)

