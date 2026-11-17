import csv

inp = "metadata/tree1217_annotation_final.tsv"
outp = "metadata/tree1217_annotation_clean.tsv"

ancestry_map = {
    "taurine": "Taurine",
    "indicine": "Indicine",
    "admixed": "Admixed",
    "review_composition": "Needs manual review",
    "metadata_unavailable": "Metadata not accessible"
}

purpose_map = {
    "beef": "Beef",
    "dairy": "Dairy",
    "dual-purpose": "Beef + dairy",
    "dairy_or_dual": "Purpose unsure",
    "crossbred": "Purpose unsure",
    "local": "Purpose unsure",
    "dairy_cross": "Purpose unsure",
    "metadata_unavailable": "Metadata not accessible"
}

composition_map = {
    "crossbred": "Crossbred",
    "dairy_cross": "Dairy cross",
    "local": "Local breed",
    "metadata_unavailable": "Metadata not accessible"
}

with open(inp) as f, open(outp, "w", newline="") as out:
    reader = csv.DictReader(f, delimiter="\t")

    fields = [
        "sample_id",
        "source",
        "genetic_ancestry",
        "production_purpose",
        "breed_composition",
        "label_flag"
    ]

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in reader:
        raw_anc = row["ancestry_group"].strip()
        raw_prod = row["production_type"].strip()

        ancestry = ancestry_map.get(
            raw_anc,
            "Needs manual review"
        )

        purpose = purpose_map.get(
            raw_prod,
            "Purpose unsure"
        )

        if raw_prod in composition_map:
            composition = composition_map[raw_prod]
        else:
            composition = "Standard breed"

        # LOCAL174 currently lacks accessible ancestry/production metadata
        if row["source"] == "LOCAL174":
            ancestry = "Metadata not accessible"
            purpose = "Metadata not accessible"
            composition = "Metadata not accessible"

        writer.writerow({
            "sample_id": row["sample_id"],
            "source": row["source"],
            "genetic_ancestry": ancestry,
            "production_purpose": purpose,
            "breed_composition": composition,
            "label_flag": row["label_flag"]
        })

print("Written:", outp)
