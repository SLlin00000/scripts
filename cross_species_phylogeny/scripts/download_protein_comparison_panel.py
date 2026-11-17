import io
import time
from pathlib import Path
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

out = Path("raw/protein_comparison")
out.mkdir(parents=True, exist_ok=True)
session = requests.Session()
records = []

ensembl = [
    ("pig_Ensembl_1865", "ENSSSCP00000062347", "sus_scrofa", 1865),
    ("sheep_main_1040", "ENSOARP00020079223", "ovis_aries", 1040),
    ("sheep_alt_1865", "ENSOARP00000004549",
     "ovis_aries_gca022416775v1", 1865),
    ("dolphin_Ensembl_591", "ENSTTRP00000016068",
     "tursiops_truncatus", 591),
]

ncbi = [
    ("cattle_NCBI_1871", "XP_059743719.1", 1871),
    ("pig_NCBI_603", "NP_001182274.1", 603),
    ("dolphin_NCBI_742", "XP_073661195.1", 742),
    ("dolphin_NCBI_603", "XP_019783990.1", 603),
    ("sheep_NCBI_1863_flagged", "XP_042107198.1", 1863),
]

def request(url, params):
    for attempt in range(3):
        try:
            time.sleep(0.5)
            r = session.get(url, params=params, timeout=(15, 60))
            r.raise_for_status()
            return r
        except requests.RequestException:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))

def save(label, accession, seq, expected):
    seq = seq.upper()
    if len(seq) != expected:
        raise ValueError(
            f"{label}: expected {expected} aa, returned {len(seq)}"
        )
    if not seq or set(seq) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO*"):
        raise ValueError(f"{label}: invalid protein sequence")
    record = SeqRecord(
        Seq(seq), id=label, description=f"accession={accession}"
    )
    SeqIO.write(record, out / f"{label}.faa", "fasta")
    records.append(record)
    print(f"{label}: {len(seq)} aa | {accession}", flush=True)

failed = []

for label, accession, species, expected in ensembl:
    try:
        r = request(
            f"https://rest.ensembl.org/sequence/id/{accession}",
            {
                "type": "protein",
                "species": species,
                "content-type": "application/json",
            }
        )
        data = r.json()
        if "error" in data:
            raise ValueError(str(data["error"]))
        (out / f"{label}.response.json").write_text(r.text)
        save(label, accession, data["seq"], expected)
    except Exception as error:
        failed.append(label)
        print(f"FAILED {label}: {error}", flush=True)

for label, accession, expected in ncbi:
    try:
        r = request(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            {
                "db": "protein", "id": accession,
                "rettype": "fasta", "retmode": "text",
                "tool": "LCORL_cross_species",
            }
        )
        record = SeqIO.read(io.StringIO(r.text), "fasta")
        (out / f"{label}.original.faa").write_text(r.text)
        save(label, accession, str(record.seq), expected)
    except Exception as error:
        failed.append(label)
        print(f"FAILED {label}: {error}", flush=True)

if failed:
    raise SystemExit(
        "Panel incomplete; failed entries: " + ", ".join(failed)
    )

panel = out / "LCORL_protein_comparison_9seq.faa"
SeqIO.write(records, panel, "fasta")
print("\nComplete comparison panel:", panel.resolve())
