import json
import time
from pathlib import Path
import requests

root = Path(__file__).resolve().parents[1]
metadata = root / "metadata" / "hybrids"
raw = root / "raw" / "hybrids"
metadata.mkdir(parents=True, exist_ok=True)
raw.mkdir(parents=True, exist_ok=True)

targets = [
    ("hybridangus", "bos_taurus_hybrid", "UOA_Angus_1"),
    ("hybridbrahman", "bos_indicus_hybrid", "UOA_Brahman_1"),
]

session = requests.Session()
session.headers.update({"Accept": "application/json"})
base = "https://rest.ensembl.org"

def get_json(endpoint, **params):
    params["content-type"] = "application/json"
    for attempt in range(3):
        try:
            time.sleep(0.5)
            response = session.get(
                base + endpoint, params=params, timeout=(15, 60)
            )
            response.raise_for_status()
            data = response.json()
            if "error" in data:
                raise RuntimeError(str(data["error"]))
            return data
        except (requests.RequestException, ValueError, RuntimeError):
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))

failed = []

for label, species, expected_assembly in targets:
    print(f"\n--- {label} ---", flush=True)
    try:
        gene = get_json(
            f"/lookup/symbol/{species}/LCORL", expand=1
        )

        if gene.get("assembly_name") != expected_assembly:
            raise RuntimeError(
                f"Assembly mismatch: {gene.get('assembly_name')}"
            )
        if gene.get("object_type") != "Gene":
            raise RuntimeError("Returned record is not a gene")

        (metadata / f"{label}.LCORL.lookup.json").write_text(
            json.dumps(gene, indent=2)
        )

        gene_id = gene["id"]
        expected_length = gene["end"] - gene["start"] + 1

        print("Assembly:", gene["assembly_name"])
        print("Gene:", gene_id, gene.get("display_name"))
        print(
            "Location (1-based, inclusive):",
            f"{gene['seq_region_name']}:{gene['start']}-{gene['end']}"
        )
        print("Strand:", gene["strand"])
        print("Expected gene length:", expected_length, "bp")

        sequence = get_json(
            f"/sequence/id/{gene_id}", type="genomic"
        )
        seq = sequence["seq"].upper()
        if len(seq) != expected_length:
            raise RuntimeError(
                f"Length mismatch: {len(seq)} != {expected_length}"
            )
        if set(seq) - set("ACGTRYSWKMBDHVN"):
            raise RuntimeError("Unexpected characters in gene sequence")

        fasta = raw / f"{label}.gene.fasta"
        with fasta.open("w") as handle:
            handle.write(
                f">{label} gene={gene_id} assembly={expected_assembly}\n"
            )
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")

        (metadata / f"{label}.gene.sequence.json").write_text(
            json.dumps(sequence, indent=2)
        )
        print("Gene FASTA saved:", fasta)
        print("N bases:", seq.count("N"))
        print("Canonical transcript:", gene.get("canonical_transcript"))

        transcripts = gene.get("Transcript", [])
        print("Transcript candidates:", len(transcripts))
        for transcript in transcripts:
            protein = transcript.get("Translation") or {}
            print(
                "Transcript:", transcript["id"],
                "| biotype:", transcript.get("biotype"),
                "| canonical:", transcript.get("is_canonical", 0),
                "| protein:", protein.get("id", "none"),
                "| AA length:", protein.get("length", "unknown"),
                "| exons:", len(transcript.get("Exon", []))
            )

    except Exception as error:
        failed.append(label)
        print(
            f"FAILED: {type(error).__name__}: {error}",
            flush=True
        )

print("\nFailed entries:", ", ".join(failed) if failed else "None")
