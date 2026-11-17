import csv
import io
import json
import time
from pathlib import Path

import requests
from Bio import SeqIO

root = Path(__file__).resolve().parents[1]
source = root / "metadata" / "ncbi_gene_candidates"
out = root / "raw" / "ncbi_genes"
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

session = requests.Session()
endpoint = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
)

def download(params):
    for attempt in range(3):
        try:
            time.sleep(0.4)
            response = session.get(
                endpoint, params=params, timeout=(15, 90)
            )
            response.raise_for_status()
            record = SeqIO.read(io.StringIO(response.text), "fasta")
            return response.text, record
        except (requests.RequestException, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))

rows = []
failed = []

for label, gene_id in targets.items():
    print(f"\n--- {label} ---", flush=True)
    try:
        data = json.loads(
            (source / f"{label}.summary.json").read_text()
        )
        gene = data["result"][gene_id]

        regions = [
            region for region in gene.get("genomicinfo", [])
            if region.get("chraccver", "").startswith("NC_")
        ]
        if len(regions) != 1:
            raise RuntimeError(
                f"Expected one chromosome location; found {len(regions)}"
            )

        region = regions[0]
        accession = region["chraccver"]

        # NCBI Gene summary coordinates are 0-based, inclusive.
        a = int(region["chrstart"])
        b = int(region["chrstop"])
        start = min(a, b) + 1
        end = max(a, b) + 1
        strand = 1 if a <= b else 2
        expected = end - start + 1

        fasta_path = out / f"{label}.gene.fasta"
        params = {
            "db": "nuccore",
            "id": accession,
            "seq_start": start,
            "seq_stop": end,
            "strand": strand,
            "rettype": "fasta",
            "retmode": "text",
            "tool": "LCORL_cross_species",
        }

        if fasta_path.exists():
            record = SeqIO.read(fasta_path, "fasta")
            provenance_path = out / f"{label}.download.json"
            if not provenance_path.exists():
                raise RuntimeError("Existing FASTA lacks download metadata")
            previous = json.loads(provenance_path.read_text())
            if previous != params:
                raise RuntimeError("Existing FASTA uses different parameters")
        else:
            original, record = download(params)
            if len(record.seq) != expected:
                raise RuntimeError(
                    f"Length mismatch: {len(record.seq)} != {expected}"
                )

            (out / f"{label}.ncbi_original.fasta").write_text(original)
            record.id = label
            record.name = label
            record.description = (
                f"gene=LCORL GeneID={gene_id} "
                f"region={accession}:{start}-{end} "
                f"strand={'+' if strand == 1 else '-'}"
            )
            SeqIO.write(record, fasta_path, "fasta")
            (out / f"{label}.download.json").write_text(
                json.dumps(params, indent=2)
            )

        seq = str(record.seq).upper()
        if len(seq) != expected:
            raise RuntimeError("Saved FASTA length mismatch")
        if set(seq) - set("ACGTRYSWKMBDHVN"):
            raise RuntimeError("Unexpected sequence characters")

        n_count = seq.count("N")
        rows.append({
            "label": label,
            "gene_id": gene_id,
            "accession": accession,
            "start_1based": start,
            "end_1based": end,
            "strand": "+" if strand == 1 else "-",
            "length_bp": len(seq),
            "N_bases": n_count,
            "N_percent": round(100 * n_count / len(seq), 4),
            "fasta": str(fasta_path.relative_to(root)),
        })
        print(
            f"{accession}:{start}-{end} | "
            f"{len(seq):,} bp | N={n_count}",
            flush=True
        )

    except Exception as error:
        failed.append(label)
        print(f"FAILED: {type(error).__name__}: {error}", flush=True)

if rows:
    report = root / "metadata" / "ncbi_gene_download_QC.tsv"
    with report.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    print("\nQC table:", report)

print("NCBI gene sequences validated:", len(rows), "/ 10")
print("Hybrid gene sequences were downloaded in the previous step.")
print("Failed entries:", ", ".join(failed) if failed else "None")
