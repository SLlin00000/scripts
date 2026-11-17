import csv
from pathlib import Path
from collections import Counter

ROOT = Path(
    "/QRISdata/Q9774/lcorl_project/"
    "LCORL_full_bam_analysis/methodology_results_reconstruction"
)
BASE = ROOT / "results"
OUT = ROOT / "presentation_tables_v2/results"

def load(path):
    assert path.is_file(), f"Missing input: {path}"
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows, f"Empty input: {path}"
    return rows

def index_unique(rows, key):
    result = {}
    for row in rows:
        value = row[key].strip()
        assert value, f"Empty {key}"
        assert value not in result, f"Duplicate {key}: {value}"
        result[value] = row
    return result

def save(path, rows, fields):
    assert not path.exists(), f"Refusing to overwrite: {path}"
    with path.open("x", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fields, delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

MISSING = {
    "", ".", "na", "n/a", "nan", "none",
    "null", "unknown", "unavailable",
    "metadata_unavailable", "not_available",
    "not_assessed"
}

def available(value):
    return value.strip().lower() not in MISSING

# ==========================================
# TABLE 1: COMPLETE METADATA FROM 1217 TREE
# ==========================================

metadata = load(
    Path("/QRISdata/Q9774/lcorl_project/"
         "cattle_phylogeny/metadata/"
         "tree1217_annotation_final.tsv")
)

assert len(metadata) == 1217
index_unique(metadata, "sample_id")

required = [
    "sample_id", "source",
    "ancestry_group", "production_type",
    "breed", "annotation_status"
]

assert all(
    key in metadata[0] for key in required
), "Metadata columns missing"

complete = []
incomplete = []

for row in metadata:
    missing = [
        key for key in required
        if not available(row[key])
    ]

    # Explicitly reject source-only records.
    if row["annotation_status"].strip().lower() == "source_only":
        missing.append("source_only_status")

    if missing:
        incomplete.append((row["sample_id"], missing))
    else:
        complete.append(row)

assert len(complete) + len(incomplete) == 1217

table1 = sorted(
    complete, key=lambda r: r["sample_id"]
)

# ==========================================
# TABLE 2: H1-H26 CORRESPONDENCE
# ==========================================

haplotypes = load(
    BASE / "tree_concordance_30018309/"
    "H1_H26_phylogenetic_correspondence.tsv"
)

branches = load(
    BASE / "branch_mapping_30019348/"
    "H1_H26_branch_correspondence_summary.tsv"
)

assignments = load(
    BASE / "run_30012993/"
    "cattle_diplotype_assignments.tsv"
)

hap_by_id = index_unique(haplotypes, "haplotype_id")
branch_by_id = index_unique(branches, "haplotype_id")
assignment_by_id = index_unique(assignments, "sample_id")

assert len(hap_by_id) == 26
assert len(branch_by_id) == 26
assert len(assignment_by_id) == 165
assert set(hap_by_id) == set(branch_by_id)
assert set(hap_by_id) == {
    f"H{i}" for i in range(1, 27)
}

copy_counts = Counter()
carrier_counts = Counter()
carrier_ids = {h: set() for h in hap_by_id}

for sample_id, row in assignment_by_id.items():
    a = row["haplotype_copy0"]
    b = row["haplotype_copy1"]

    assert a in hap_by_id and b in hap_by_id

    copy_counts.update([a, b])
    carrier_counts.update({a, b})

    carrier_ids[a].add(sample_id)
    carrier_ids[b].add(sample_id)

assert sum(copy_counts.values()) == 330

table2 = []

for i in range(1, 27):
    hid = f"H{i}"
    h = hap_by_id[hid]
    b = branch_by_id[hid]

    assert (
        h["amino_acid_haplotype"]
        == b["amino_acid_haplotype"]
    )
    assert int(h["chromosome_count"]) == copy_counts[hid]
    assert int(h["carrier_count"]) == carrier_counts[hid]
    assert int(b["chromosome_count"]) == copy_counts[hid]
    assert int(b["carrier_count"]) == carrier_counts[hid]

    assert abs(
        float(h["chromosome_frequency"])
        - copy_counts[hid] / 330
    ) < 1e-6

    assert abs(
        float(h["carrier_frequency"])
        - carrier_counts[hid] / 165
    ) < 1e-6

    source_tips = set(
        x for x in h["tree_tip_ids"].split(";") if x
    )
    assert source_tips == carrier_ids[hid], (
        f"Tree-tip mismatch for {hid}"
    )

    table2.append({
        "haplotype_id": hid,
        "amino_acid_haplotype":
            h["amino_acid_haplotype"],
        "ALT_missense_changes":
            h["ALT_missense_changes"],
        "chromosome_count":
            h["chromosome_count"],
        "chromosome_frequency":
            h["chromosome_frequency"],
        "carrier_count":
            h["carrier_count"],
        "carrier_frequency":
            h["carrier_frequency"],
        "diplotype_combinations":
            h["diplotype_combinations"],
        "tree_tip_ids":
            h["tree_tip_ids"],
        "selected_split_id":
            b["selected_split_id"],
        "selected_side":
            b["selected_side"],
        "assessed_on_selected_side":
            b["assessed_on_selected_side"],
        "carriers_on_selected_side":
            b["carriers_on_selected_side"],
        "carrier_proportion_on_selected_side":
            b["carrier_proportion_on_selected_side"],
        "carrier_coverage_on_selected_side":
            b["carrier_coverage_on_selected_side"],
        "proportion_difference":
            b["proportion_difference"],
        "selection_status":
            b["selection_status"]
    })

assert len(table2) == 26

# ==========================================
# WRITE NEW TABLES ONLY AFTER VALIDATION
# ==========================================

OUT.mkdir(parents=True, exist_ok=True)

t1 = OUT / "Table1_complete_metadata_tree1217.tsv"
t2 = OUT / "Table2_H1_H26_tree_correspondence.tsv"

assert not t1.exists()
assert not t2.exists()

save(t1, table1, list(metadata[0].keys()))
save(t2, table2, list(table2[0].keys()))

print("===== TABLE 1: COMPLETE METADATA =====")
print("Original tree samples:", len(metadata))
print("Complete metadata samples:", len(table1))
print("Excluded incomplete samples:", len(incomplete))
print("Included sources:", dict(Counter(
    r["source"] for r in table1
)))
print("Excluded sources:", dict(Counter(
    r["source"] for r in metadata
    if r["sample_id"] in {
        sid for sid, _ in incomplete
    }
)))

print("\nMissing/incomplete reasons:")
for reason, n in sorted(
    Counter(
        reason
        for _, reasons in incomplete
        for reason in set(reasons)
    ).items()
):
    print(f"{reason}: {n}")

print("\n===== TABLE 2: HAPLOTYPES =====")
print("Protein haplotypes:", len(table2))
print("Phased cattle:", len(assignments))
print("Chromosome copies:", sum(copy_counts.values()))
print("H1 copies:", copy_counts["H1"])
print("H1 carriers:", carrier_counts["H1"])

print("\n===== OUTPUT FILES =====")
print(t1)
print(t2)
print("VALIDATION PASSED")
