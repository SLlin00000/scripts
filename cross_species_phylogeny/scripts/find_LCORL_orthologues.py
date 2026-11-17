import json
import time
from pathlib import Path
import requests

out = Path("metadata/source_audit/ensembl_comparison")
out.mkdir(parents=True, exist_ok=True)

base = "https://rest.ensembl.org"
session = requests.Session()
source_id = "ENSSSCG00000008748"

def query(endpoint, **params):
    params["content-type"] = "application/json"
    for attempt in range(3):
        time.sleep(0.5)
        try:
            response = session.get(
                base + endpoint, params=params, timeout=(15, 90)
            )
            if response.status_code == 400:
                raise RuntimeError(response.text[:1000])
            response.raise_for_status()
            data = response.json()
            if isinstance(data, dict) and "error" in data:
                raise RuntimeError(str(data["error"]))
            return data
        except (requests.RequestException, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))

for species in ["ovis_aries", "tursiops_truncatus"]:
    print(f"\n--- {species} ---", flush=True)
    try:
        data = query(
            f"/homology/id/sus_scrofa/{source_id}",
            target_species=species,
            type="orthologues"
        )
        (out / f"{species}.pig_LCORL_orthologues.json").write_text(
            json.dumps(data, indent=2)
        )

        hits = [
            h for entry in data.get("data", [])
            for h in entry.get("homologies", [])
        ]
        print("Orthology records:", len(hits))

        seen = set()
        for hit in hits:
            target = hit["target"]
            gene_id = target["id"]
            print(
                "Target gene:", gene_id,
                "| Species:", target.get("species"),
                "| Orthology:", hit.get("type"),
                "| Confidence:", hit.get("is_high_confidence", "not returned")
            )

            if target.get("species") != species or gene_id in seen:
                continue
            seen.add(gene_id)

            gene = query(
                f"/lookup/id/{gene_id}", species=species, expand=1
            )
            (out / f"{species}.{gene_id}.lookup.json").write_text(
                json.dumps(gene, indent=2)
            )

            print("Symbol:", gene.get("display_name"))
            print("Description:", gene.get("description"))
            print("Assembly:", gene.get("assembly_name"))
            print(
                "Location:",
                gene.get("seq_region_name"),
                gene.get("start"), gene.get("end"),
                "| Strand:", gene.get("strand")
            )

            for transcript in gene.get("Transcript", []):
                protein = transcript.get("Translation") or {}
                if protein.get("id"):
                    print(
                        f"Transcript={transcript['id']} | "
                        f"Protein={protein['id']} | "
                        f"AA={protein.get('length')} | "
                        f"Exons={len(transcript.get('Exon', []))} | "
                        f"Canonical={transcript.get('is_canonical', 0)}"
                    )

    except Exception as error:
        print(f"FAILED: {type(error).__name__}: {error}", flush=True)
