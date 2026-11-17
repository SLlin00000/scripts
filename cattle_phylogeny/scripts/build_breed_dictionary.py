#!/usr/bin/env python3

import csv

INFILE = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/run8_unique_breed_candidates.txt"
OUTFILE = "/QRISdata/Q9774/lcorl_project/cattle_phylogeny/metadata/breed_dictionary.tsv"


# raw label :
# (normalised breed, ancestry expectation, production type, status, note)

MAPPING = {

    # -------------------------
    # DAIRY
    # -------------------------
    "Holstein":
        ("Holstein", "taurine", "dairy", "KEEP", ""),

    "Holstein bull":
        ("Holstein", "taurine", "dairy", "KEEP", "breed label normalised"),

    "Polish Holstein-Friesian":
        ("Holstein-Friesian", "taurine", "dairy", "KEEP", ""),

    "Jersey":
        ("Jersey", "taurine", "dairy", "KEEP", ""),

    "Norwegian Red":
        ("Norwegian Red", "taurine", "dairy", "KEEP", ""),

    "Friesian-Jersey cross":
        ("Friesian x Jersey", "taurine", "dairy_cross", "KEEP", "crossbred dairy"),


    # -------------------------
    # BEEF TAURINE
    # -------------------------
    "Angus":
        ("Angus", "taurine", "beef", "KEEP", ""),

    "German Angus":
        ("German Angus", "taurine", "beef", "KEEP", ""),

    "Red Angus":
        ("Red Angus", "taurine", "beef", "KEEP", ""),

    "Angus lowline":
        ("Lowline Angus", "taurine", "beef", "KEEP", "normalised name"),

    "Hereford":
        ("Hereford", "taurine", "beef", "KEEP", ""),

    "Herford":
        ("Hereford", "taurine", "beef", "KEEP", "spelling normalised"),

    "Hereford miniature":
        ("Miniature Hereford", "taurine", "beef", "KEEP", ""),

    "Charolais":
        ("Charolais", "taurine", "beef", "KEEP", ""),

    "Limousin":
        ("Limousin", "taurine", "beef", "KEEP", ""),

    "Belgian Blue":
        ("Belgian Blue", "taurine", "beef", "KEEP", ""),

    "Salers":
        ("Salers", "taurine", "beef", "KEEP", ""),

    "Piedmontese":
        ("Piedmontese", "taurine", "beef", "KEEP", ""),

    "Chianina":
        ("Chianina", "taurine", "beef", "KEEP", ""),

    "Texas Longhorn":
        ("Texas Longhorn", "taurine", "beef", "KEEP", ""),

    "Belted Galloway":
        ("Belted Galloway", "taurine", "beef", "KEEP", ""),

    "Highland":
        ("Highland", "taurine", "beef", "KEEP", ""),

    "highland":
        ("Highland", "taurine", "beef", "KEEP", "capitalisation normalised"),

    "Hanwoo":
        ("Hanwoo", "taurine", "beef", "KEEP", ""),

    "Corriente":
        ("Corriente", "taurine", "beef", "KEEP", ""),

    "Maine Anjou":
        ("Maine-Anjou", "taurine", "beef", "KEEP", "name normalised"),

    "Maine-Anjou":
        ("Maine-Anjou", "taurine", "beef", "KEEP", ""),

    "Blonde d'Aquitaine":
        ("Blonde d'Aquitaine", "taurine", "beef", "KEEP", ""),

    "Dexter":
        ("Dexter", "taurine", "beef", "KEEP", ""),


    # -------------------------
    # DUAL PURPOSE
    # -------------------------
    "Simmental":
        ("Simmental", "taurine", "dual-purpose", "KEEP", ""),

    "Fleckvieh":
        ("Fleckvieh", "taurine", "dual-purpose", "KEEP", ""),

    "Gelbvieh":
        ("Gelbvieh", "taurine", "dual-purpose", "KEEP", ""),

    "Braunvieh":
        ("Braunvieh", "taurine", "dual-purpose", "KEEP", ""),

    "Braunveih":
        ("Braunvieh", "taurine", "dual-purpose", "KEEP", "spelling normalised"),

    "Shorthorn":
        ("Shorthorn", "taurine", "dual-purpose", "KEEP", ""),


    # -------------------------
    # INDICINE / ZEBU
    # -------------------------
    "Brahman":
        ("Brahman", "indicine", "beef", "KEEP", ""),

    "Nelore":
        ("Nelore", "indicine", "beef", "KEEP", ""),

    "Gir":
        ("Gir", "indicine", "dairy_or_dual", "KEEP", ""),

    "Tharparkar":
        ("Tharparkar", "indicine", "dual-purpose", "KEEP", ""),

    "Indian Zebu":
        ("Indian Zebu", "indicine", "unknown", "REVIEW", "broad population label"),

    "East African Shorthorn Zebu":
        ("East African Shorthorn Zebu", "indicine", "local", "KEEP", ""),

    "Small East African Zebu":
        ("Small East African Zebu", "indicine", "local", "KEEP", ""),

    "Shorthorn Zebu":
        ("Shorthorn Zebu", "indicine", "local", "KEEP", ""),

    "Mangshi Zebu":
        ("Mangshi Zebu", "indicine", "local", "KEEP", ""),


    # -------------------------
    # KNOWN COMPOSITE / CROSS
    # -------------------------
    "Brangus":
        ("Brangus", "admixed", "beef", "KEEP", ""),

    "Beefmaster":
        ("Beefmaster", "admixed", "beef", "KEEP", ""),

    "Santa Gertrudis":
        ("Santa Gertrudis", "admixed", "beef", "KEEP", ""),

    "Ankole Zebu cross":
        ("Ankole Zebu cross", "admixed", "crossbred", "KEEP", ""),

    "Holstein x Hereford":
        ("Holstein x Hereford", "taurine", "crossbred", "KEEP", ""),

    "Limousin x Hereford":
        ("Limousin x Hereford", "taurine", "crossbred", "KEEP", ""),

    "Charolais x Angus":
        ("Charolais x Angus", "taurine", "crossbred", "KEEP", ""),

    "Charolais-Red Angus":
        ("Charolais x Red Angus", "taurine", "crossbred", "KEEP", ""),

    "BelgianBlue/LimousinCross":
        ("Belgian Blue x Limousin", "taurine", "crossbred", "KEEP", ""),


    # -------------------------
    # NOT RELIABLE BREED LABELS
    # -------------------------
    "sleeping beauty transposon":
        ("", "unknown", "unknown", "EXCLUDE", "experimental construct label"),

    "piggybac transposon":
        ("", "unknown", "unknown", "EXCLUDE", "experimental construct label"),

    "Cloned bull":
        ("", "unknown", "unknown", "EXCLUDE", "not a breed"),

    "F1":
        ("", "unknown", "crossbred", "REVIEW", "parental breeds unspecified"),

    "Composite":
        ("", "admixed", "unknown", "REVIEW", "composition unspecified"),

    "CPC98":
        ("", "unknown", "unknown", "REVIEW", "coded label"),

    "DRS001209":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001210":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001211":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001212":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001213":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001214":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001215":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001216":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001217":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001218":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),

    "DRS001219":
        ("", "unknown", "unknown", "EXCLUDE", "accession-like label"),
}


with open(INFILE) as inp, open(OUTFILE, "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "raw_breed",
        "normalised_breed",
        "expected_ancestry",
        "production_type",
        "status",
        "note"
    ])

    for line in inp:
        raw = line.strip()

        if not raw:
            continue

        if raw in MAPPING:
            norm, anc, prod, status, note = MAPPING[raw]

        elif raw.startswith("Crossbred"):
            norm = raw
            anc = "review_composition"
            prod = "crossbred"
            status = "KEEP"
            note = "explicit crossbred composition; ancestry to be derived from parents"

        else:
            norm = raw
            anc = ""
            prod = ""
            status = "REVIEW"
            note = "manual verification required"

        writer.writerow([
            raw,
            norm,
            anc,
            prod,
            status,
            note
        ])

