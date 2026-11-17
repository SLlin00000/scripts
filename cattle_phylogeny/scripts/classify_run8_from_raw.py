#!/usr/bin/env python3

import csv
import re

INFILE = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/run8_metadata_raw.tsv"
OUTFILE = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/run8_metadata_reclassified.tsv"

MISSING = {
    "", "na", "n/a", "none", "null", "missing",
    "unknown", "not provided", "not collected"
}

def clean(x):
    return (x or "").strip()

def valid(x):
    return clean(x).lower() not in MISSING

def ancestry_from_scientific_name(name):
    s = clean(name).lower()

    if "bos indicus" in s and "bos taurus" in s:
        return "admixed"
    if "bos indicus" in s:
        return "indicine"
    if "bos taurus" in s:
        return "taurine"
    if "bos primigenius" in s:
        return "review"

    return "unknown"

# Conservative known breed names/aliases for rescue from description
ALIASES = {
    "holstein": "Holstein",
    "brown swiss": "Brown Swiss",
    "brown suisse": "Brown Swiss",
    "braunvieh": "Braunvieh",
    "original braunvieh": "Braunvieh",
    "simmental": "Simmental",
    "angus": "Angus",
    "red angus": "Red Angus",
    "hereford": "Hereford",
    "charolais": "Charolais",
    "limousin": "Limousin",
    "jersey": "Jersey",
    "fleckvieh": "Fleckvieh",
    "gelbvieh": "Gelbvieh",
    "belgian blue": "Belgian Blue",
    "tyrolean grey cattle": "Tyrolean Grey",
    "tyrolean grey": "Tyrolean Grey",
    "brahman": "Brahman",
    "nelore": "Nelore",
}

ALIASES_SORTED = sorted(ALIASES.keys(), key=len, reverse=True)

def rescue_from_description(description):
    d = clean(description).lower()

    if not d:
        return "", "", "unresolved"

    matches = []

    for alias in ALIASES_SORTED:
        pattern = r"(?<![A-Za-z])" + re.escape(alias) + r"(?![A-Za-z])"

        if re.search(pattern, d, flags=re.IGNORECASE):
            norm = ALIASES[alias]
            if norm not in matches:
                matches.append(norm)

    if len(matches) == 1:
        return matches[0], "description", "medium"

    if len(matches) > 1:
        return " x ".join(matches), "description_cross", "medium"

    return "", "", "unresolved"


with open(INFILE, newline="") as inp, open(OUTFILE, "w", newline="") as out:

    reader = csv.DictReader(inp, delimiter="\t")

    fields = reader.fieldnames + [
        "breed_candidate",
        "breed_source",
        "breed_confidence",
        "ancestry_group",
        "label_status"
    ]

    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for row in reader:

        ancestry = ancestry_from_scientific_name(row["scientific_name"])

        breed_candidate = ""
        breed_source = ""
        confidence = "unresolved"

        # Priority 1: explicit breed
        if valid(row["breed"]):
            breed_candidate = clean(row["breed"])
            breed_source = "breed"
            confidence = "high"

        # Priority 2: informative isolate
        elif valid(row["isolate"]):
            breed_candidate = clean(row["isolate"])
            breed_source = "isolate"
            confidence = "review"

        # Priority 3: description rescue
        else:
            rescued, source, conf = rescue_from_description(row["description"])

            if rescued:
                breed_candidate = rescued
                breed_source = source
                confidence = conf

        if confidence == "high":
            status = "well_labelled"

        elif confidence in {"medium", "review"}:
            status = "needs_review"

        else:
            status = "poorly_labelled"

        row["breed_candidate"] = breed_candidate
        row["breed_source"] = breed_source
        row["breed_confidence"] = confidence
        row["ancestry_group"] = ancestry
        row["label_status"] = status

        writer.writerow(row)

