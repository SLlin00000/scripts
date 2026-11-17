import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path
import requests

root = Path(__file__).resolve().parents[1]
out = root / "raw" / "ncbi_records"
out.mkdir(parents=True, exist_ok=True)

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

base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
session = requests.Session()

def fetch(endpoint, **params):
    params["tool"] = "LCORL_cross_species"
    for attempt in range(3):
        try:
            time.sleep(0.4)
            r = session.get(
                base + endpoint, params=params, timeout=(15, 60)
            )
            r.raise_for_status()
            return r.text
        except requests.RequestException:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))

failed = []

for label, gene_id in targets.items():
    print(f"\n--- {label} | GeneID {gene_id} ---", flush=True)
    try:
        path = out / f"{label}.gene.xml"
        if path.exists():
            text = path.read_text()
        else:
            text = fetch(
                "efetch.fcgi", db="gene", id=gene_id, retmode="xml"
            )

        tree = ET.fromstring(text)
        if tree.find(".//ERROR") is not None:
            raise RuntimeError(tree.findtext(".//ERROR"))
        if tree.find(".//Gene-track_geneid") is None:
            raise RuntimeError("No gene record returned")

        path.write_text(text)

        headings = sorted({
            x.text for x in tree.findall(".//Gene-commentary_heading")
            if x.text
        })
        print("Record headings:", "; ".join(headings))

        seen = set()
        for item in tree.iter("Gene-commentary"):
            typ = item.find("Gene-commentary_type")
            kind = typ.get("value", "") if typ is not None else ""
            if kind not in {"mRNA", "peptide"}:
                continue

            acc = item.findtext("Gene-commentary_accession", "")
            version = item.findtext("Gene-commentary_version", "")
            accession = f"{acc}.{version}" if version else acc
            label_text = item.findtext("Gene-commentary_label", "")
            product = item.findtext("Gene-commentary_text", "")
            row = (kind, accession, label_text, product)

            if acc and row not in seen:
                seen.add(row)
                print(" | ".join(row))

        if not seen:
            print("No mRNA/protein accessions found in this record.")

    except Exception as error:
        failed.append(label)
        print(f"FAILED: {type(error).__name__}: {error}", flush=True)

for assembly in ["UOA_Angus_1", "UOA_Brahman_1"]:
    print(f"\n--- Assembly: {assembly} ---", flush=True)
    try:
        search_text = fetch(
            "esearch.fcgi", db="assembly",
            term=f'"{assembly}"[Assembly Name]',
            retmode="json", retmax=10
        )
        search = json.loads(search_text)
        (out / f"{assembly}.search.json").write_text(search_text)
        if "error" in search:
            raise RuntimeError(str(search["error"]))

        result = search["esearchresult"]
        print("Matching assemblies:", result["count"])
        ids = result["idlist"]
        if not ids:
            continue

        summary_text = fetch(
            "esummary.fcgi", db="assembly",
            id=",".join(ids), retmode="json"
        )
        summary = json.loads(summary_text)
        (out / f"{assembly}.summary.json").write_text(summary_text)
        if "error" in summary:
            raise RuntimeError(str(summary["error"]))

        for uid in summary["result"]["uids"]:
            record = summary["result"][uid]
            for key in [
                "assemblyname", "assemblyaccession", "organism",
                "assemblystatus", "synonym",
                "ftppath_refseq", "ftppath_genbank"
            ]:
                print(f"{key}: {record.get(key, '')}")

    except Exception as error:
        failed.append(assembly)
        print(f"FAILED: {type(error).__name__}: {error}", flush=True)

print("\nSaved records:", out)
print("Failed queries:", ", ".join(failed) if failed else "None")
