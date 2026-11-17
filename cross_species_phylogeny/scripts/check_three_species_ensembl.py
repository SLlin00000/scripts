import json
import time
from pathlib import Path
import requests

root = Path(__file__).resolve().parents[1]
out = root / "metadata/source_audit/ensembl_comparison"
out.mkdir(parents=True, exist_ok=True)

targets = [
    ("pig", "sus_scrofa"),
    ("dolphin", "tursiops_truncatus"),
    ("sheep", "ovis_aries"),
]

session = requests.Session()
base = "https://rest.ensembl.org"

def get_json(endpoint, **params):
    params["content-type"] = "application/json"
    for attempt in range(3):
        try:
            time.sleep(0.5)
            r = session.get(
                base + endpoint, params=params, timeout=(15, 60)
            )
            r.raise_for_status()
            data = r.json()
            if "error" in data:
                raise RuntimeError(str(data["error"]))
            return data
        except (requests.RequestException, ValueError, RuntimeError):
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))

for label, species in targets:
    print(f"\n--- {label} ---", flush=True)
    try:
        gene = get_json(
            f"/lookup/symbol/{species}/LCORL", expand=1
        )
        assembly = get_json(f"/info/assembly/{species}")

        (out / f"{label}.lookup.json").write_text(
            json.dumps(gene, indent=2)
        )
        (out / f"{label}.assembly.json").write_text(
            json.dumps(assembly, indent=2)
        )

        print("Gene:", gene["id"], gene.get("display_name"))
        print("Assembly:", gene.get("assembly_name"))
        print("Assembly accession:", assembly.get("assembly_accession"))
        print("Canonical:", gene.get("canonical_transcript"))

        count = 0
        for transcript in gene.get("Transcript", []):
            protein = transcript.get("Translation") or {}
            if not protein.get("id"):
                continue
            count += 1
            print(
                f"Transcript={transcript['id']} | "
                f"Protein={protein['id']} | "
                f"AA={protein.get('length', 'unknown')} | "
                f"Exons={len(transcript.get('Exon', []))} | "
                f"Canonical={transcript.get('is_canonical', 0)}"
            )
        print("Protein-coding candidates:", count)

    except Exception as error:
        print(f"FAILED: {type(error).__name__}: {error}", flush=True)

print("\nSaved:", out)
