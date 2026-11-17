import csv
import io
import json
import time
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition

root = Path(__file__).resolve().parents[1]
audit = root / "metadata" / "source_audit"
audit.mkdir(parents=True, exist_ok=True)

targets = {
    "pig": "100337668",
    "alpaca": "102534421",
    "camel": "105098353",
    "dolphin": "101323726",
    "sheep": "101104320",
    "goat": "102190595",
    "bostaurus": "540095",
    "bosindicus": "109560356",
    "buffalo": "102395331",
    "horse": "100068577",
}

session = requests.Session()
endpoint = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
summary = []
cds_rows = []

def fetch_record(params):
    for attempt in range(3):
        try:
            time.sleep(0.4)
            response = session.get(
                endpoint, params=params, timeout=(15, 90)
            )
            response.raise_for_status()
            record = SeqIO.read(
                io.StringIO(response.text), "genbank"
            )
            return response.text, record
        except (requests.RequestException, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))

def matches_lcorl(feature, gene_id):
    names = [x.upper() for x in feature.qualifiers.get("gene", [])]
    refs = feature.qualifiers.get("db_xref", [])
    return "LCORL" in names or f"GeneID:{gene_id}" in refs

for label, gene_id in targets.items():
    print(f"\n--- {label} ---", flush=True)
    try:
        folder = root / "raw" / "ncbi_genes"
        params = json.loads(
            (folder / f"{label}.download.json").read_text()
        )
        params.update(rettype="gbwithparts", retmode="text")

        gb_path = audit / f"{label}.region.gb"
        if gb_path.exists():
            record = SeqIO.read(gb_path, "genbank")
        else:
            text, record = fetch_record(params)
            gb_path.write_text(text)

        original = SeqIO.read(
            folder / f"{label}.gene.fasta", "fasta"
        )
        sequence_match = (
            str(original.seq).upper() == str(record.seq).upper()
        )

        gene_features = [
            f for f in record.features
            if f.type == "gene" and matches_lcorl(f, gene_id)
        ]
        cds_features = [
            f for f in record.features
            if f.type == "CDS" and matches_lcorl(f, gene_id)
        ]

        lengths = []
        partial_count = 0
        for f in cds_features:
            q = f.qualifiers
            protein = q.get("protein_id", [""])[0]
            translation = q.get("translation", [""])[0]
            parts = list(f.location.parts)
            partial = (
                "partial" in q
                or any(
                    not isinstance(edge, ExactPosition)
                    for part in parts
                    for edge in (part.start, part.end)
                )
            )
            partial_count += int(partial)
            if translation:
                lengths.append(len(translation))

            cds_rows.append({
                "label": label,
                "source": "NCBI",
                "transcript_id": q.get("transcript_id", [""])[0],
                "protein_id": protein,
                "aa_length": len(translation) if translation else "",
                "partial": partial,
                "CDS_location_in_downloaded_region": str(f.location),
                "notes": json.dumps({
                    k: q[k] for k in
                    ["product", "inference", "experiment",
                     "exception", "note", "pseudo"]
                    if k in q
                }, ensure_ascii=False),
            })

        passed = sequence_match and bool(gene_features)
        status = "PASS_IDENTITY" if passed else "REVIEW"
        summary.append({
            "label": label,
            "status": status,
            "sequence_matches_database": sequence_match,
            "LCORL_gene_features": len(gene_features),
            "CDS_candidates": len(cds_features),
            "partial_CDS_candidates": partial_count,
        })

        print("Sequence matches database:", sequence_match)
        print("LCORL gene features:", len(gene_features))
        print("CDS candidates:", len(cds_features))
        print("Partial CDS candidates:", partial_count)
        print("Distinct AA lengths:", sorted(set(lengths)))
        print("Status:", status)

    except Exception as error:
        summary.append({"label": label, "status": "FAILED"})
        print(f"FAILED: {type(error).__name__}: {error}", flush=True)

for label, expected_assembly in [
    ("hybridangus", "UOA_Angus_1"),
    ("hybridbrahman", "UOA_Brahman_1"),
]:
    print(f"\n--- {label} ---", flush=True)
    try:
        gene = json.loads(
            (root / "metadata/hybrids" /
             f"{label}.LCORL.lookup.json").read_text()
        )
        sequence = json.loads(
            (root / "metadata/hybrids" /
             f"{label}.gene.sequence.json").read_text()
        )
        fasta = SeqIO.read(
            root / "raw/hybrids" / f"{label}.gene.fasta", "fasta"
        )

        sequence_match = (
            str(fasta.seq).upper() == sequence["seq"].upper()
        )
        passed = (
            sequence_match
            and gene.get("display_name", "").upper() == "LCORL"
            and gene.get("assembly_name") == expected_assembly
            and len(fasta.seq) == gene["end"] - gene["start"] + 1
        )

        candidates = []
        for transcript in gene.get("Transcript", []):
            protein = transcript.get("Translation") or {}
            if protein.get("id"):
                candidates.append(protein)
                cds_rows.append({
                    "label": label,
                    "source": "Ensembl",
                    "transcript_id": transcript["id"],
                    "protein_id": protein["id"],
                    "aa_length": protein.get("length", ""),
                    "partial": "not_checked",
                    "CDS_location_in_downloaded_region": "",
                    "notes": json.dumps({
                        "canonical": transcript.get("is_canonical", 0)
                    }),
                })

        status = "PASS_IDENTITY" if passed else "REVIEW"
        summary.append({
            "label": label,
            "status": status,
            "sequence_matches_database": sequence_match,
            "LCORL_gene_features": 1 if passed else "",
            "CDS_candidates": len(candidates),
            "partial_CDS_candidates": "not_checked",
        })
        print("Assembly:", gene.get("assembly_name"))
        print("Gene:", gene.get("display_name"), gene["id"])
        print("Sequence matches saved database response:", sequence_match)
        print("Status:", status)

    except Exception as error:
        summary.append({"label": label, "status": "FAILED"})
        print(f"FAILED: {type(error).__name__}: {error}", flush=True)

for filename, rows in [
    ("gene_source_audit.tsv", summary),
    ("CDS_protein_candidates.tsv", cds_rows),
]:
    if rows:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with (audit / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)

print("\nReports saved in:", audit)
print("PASS_IDENTITY confirms identity/extraction, not a confidence ranking.")
