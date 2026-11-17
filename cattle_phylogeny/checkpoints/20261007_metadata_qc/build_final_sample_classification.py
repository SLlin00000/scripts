#!/usr/bin/env python3

import csv
from collections import Counter

BASE = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny"

SAMPLES = f"{BASE}/metadata/run8_metadata_reclassified.tsv"
DICT = f"{BASE}/metadata/breed_dictionary_master.tsv"
OUT = f"{BASE}/metadata/run8_samples_final_classification.tsv"


breed_map = {}

with open(DICT, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        # Ignore old dictionary labels not present in current analysis
        if row["present_in_current"].strip().lower() != "yes":
            continue

        raw = row["raw_breed"].strip()

        breed_map[raw] = {
            "normalised_breed": row["normalised_breed"].strip(),
            "expected_ancestry": row["expected_ancestry"].strip(),
            "production_type": row["production_type"].strip(),
            "dictionary_status": row["status"].strip(),
            "dictionary_note": row["note"].strip()
        }


def ancestry_qc(sample_ancestry, expected_ancestry):

    s = (sample_ancestry or "").strip().lower()
    e = (expected_ancestry or "").strip().lower()

    if not e or e in {
        "unknown",
        "review",
        "review_composition"
    }:
        return "NOT_ASSESSABLE"

    if not s or s == "unknown":
        return "NOT_ASSESSABLE"

    if s == e:
        return "CONSISTENT"

    # An admixed taxonomic label can be compatible with a breed
    # classified as taurine- or indicine-derived, and vice versa.
    if "admixed" in {s, e}:
        return "COMPATIBLE_ADMIXTURE"

    # Only direct taurine-vs-indicine disagreement is treated
    # as a true ancestry conflict.
    if {s, e} == {"taurine", "indicine"}:
        return "TRUE_CONFLICT"

    return "NOT_ASSESSABLE"


def final_status(
    breed_candidate,
    sample_label_status,
    dictionary_status,
    ancestry_check
):

    # No usable breed at all
    if not breed_candidate:
        return "EXCLUDE_POOR_METADATA"

    # Explicitly bad/non-breed dictionary labels
    if dictionary_status == "EXCLUDE":
        return "EXCLUDE_DICTIONARY"

    # Breed exists but still needs manual verification
    if dictionary_status == "REVIEW":
        return "REVIEW"

    # Dictionary says KEEP, but ancestry metadata disagree
    if ancestry_check == "TRUE_CONFLICT":
        return "REVIEW_CONFLICT"

    # Dictionary says KEEP and no ancestry problem
    if dictionary_status == "KEEP":
        return "KEEP"

    return "REVIEW"


with open(SAMPLES, newline="") as inp, open(
    OUT, "w", newline=""
) as out:

    reader = csv.DictReader(inp, delimiter="\t")

    extra_fields = [
        "normalised_breed",
        "expected_ancestry",
        "production_type",
        "dictionary_status",
        "dictionary_note",
        "ancestry_qc",
        "final_status"
    ]

    writer = csv.DictWriter(
        out,
        fieldnames=reader.fieldnames + extra_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in reader:

        breed_candidate = row["breed_candidate"].strip()
        sample_ancestry = row["ancestry_group"].strip()

        if breed_candidate in breed_map:

            info = breed_map[breed_candidate]

            normalised_breed = info["normalised_breed"]
            expected_ancestry = info["expected_ancestry"]
            production_type = info["production_type"]
            dictionary_status = info["dictionary_status"]
            dictionary_note = info["dictionary_note"]

        else:

            normalised_breed = ""
            expected_ancestry = ""
            production_type = ""
            dictionary_status = (
                "NO_BREED" if not breed_candidate else "UNMAPPED"
            )
            dictionary_note = ""

        ancestry_check = ancestry_qc(
            sample_ancestry,
            expected_ancestry
        )

        status = final_status(
            breed_candidate,
            row["label_status"].strip(),
            dictionary_status,
            ancestry_check
        )

        row["normalised_breed"] = normalised_breed
        row["expected_ancestry"] = expected_ancestry
        row["production_type"] = production_type
        row["dictionary_status"] = dictionary_status
        row["dictionary_note"] = dictionary_note
        row["ancestry_qc"] = ancestry_check
        row["final_status"] = status

        writer.writerow(row)


status_counts = Counter()
ancestry_counts = Counter()
dictionary_counts = Counter()
source_counts = Counter()

with open(OUT, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    total = 0

    for row in reader:
        total += 1
        status_counts[row["final_status"]] += 1
        ancestry_counts[row["ancestry_qc"]] += 1
        dictionary_counts[row["dictionary_status"]] += 1

        source = row["breed_source"].strip()
        if not source:
            source = "none"

        source_counts[source] += 1


print("TOTAL_SAMPLES", total)

print("\nFINAL_STATUS")
for key in sorted(status_counts):
    print(key, status_counts[key])

print("\nANCESTRY_QC")
for key in sorted(ancestry_counts):
    print(key, ancestry_counts[key])

print("\nDICTIONARY_STATUS")
for key in sorted(dictionary_counts):
    print(key, dictionary_counts[key])

print("\nBREED_SOURCE")
for key in sorted(source_counts):
    print(key, source_counts[key])

