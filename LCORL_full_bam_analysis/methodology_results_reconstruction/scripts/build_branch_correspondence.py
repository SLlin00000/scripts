#!/usr/bin/env python3
import argparse
import csv
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from Bio import Phylo

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=columns, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    project = Path(args.project).resolve()
    work = project / (
        "LCORL_full_bam_analysis/"
        "methodology_results_reconstruction"
    )
    original = work / "results/run_30012993"
    previous = work / "results/tree_concordance_30018309"
    out = Path(args.out).resolve()

    require(out.is_dir() and not any(out.iterdir()),
            "Output directory must be empty")

    for p in (original, previous):
        require(
            (p / "COMPLETED_VALIDATED.txt")
            .read_text().strip() == "PASS",
            f"Previous analysis not validated: {p}"
        )

    individuals = read_tsv(
        previous / "individual_tree_mapping.tsv"
    )
    catalogue = read_tsv(
        original / "protein_haplotype_catalogue.tsv"
    )

    require(len(individuals) == 165,
            "Expected 165 mapped cattle")
    require(len(catalogue) == 26,
            "Expected 26 protein haplotypes")

    sample_ids = set()
    carriers = defaultdict(set)

    for r in individuals:
        sample = r["tree_tip_id"]
        require(sample not in sample_ids,
                "Duplicate sample ID")
        sample_ids.add(sample)
        carriers[r["haplotype_copy0"]].add(sample)
        carriers[r["haplotype_copy1"]].add(sample)

    hap_ids = {
        r["haplotype_id"] for r in catalogue
    }
    require(set(carriers) == hap_ids,
            "Haplotype IDs do not match catalogue")

    tree_path = project / (
        "cattle_phylogeny/results/"
        "cattle1217_LCORL_NJ_tree.newick"
    )
    tree = Phylo.read(str(tree_path), "newick")

    terminals = tree.get_terminals()
    tip_names = [t.name for t in terminals]

    require(
        len(tip_names) == 1217
        and len(set(tip_names)) == 1217,
        "Invalid NJ tree tips"
    )
    require(sample_ids <= set(tip_names),
            "Mapped samples missing from tree")

    print("STEP 1: Deriving unrooted branch splits",
          flush=True)

    # Terminal sets on one side of each edge.
    descendant = {}

    def collect(node):
        if node.is_terminal():
            names = frozenset([node.name])
        else:
            names = frozenset().union(
                *(collect(c) for c in node.clades)
            )
        descendant[node] = names
        return names

    all_tips = collect(tree.root)
    require(len(all_tips) == 1217,
            "Unexpected tree tip count")

    # An unrooted split is identical when side A/B swap.
    # Canonicalize using the smaller full-tree tip set.
    unique_splits = {}

    for node in tree.find_clades(order="preorder"):
        if node is tree.root:
            continue
        side = descendant[node]
        other = all_tips - side

        if min(len(side), len(other)) < 2:
            continue

        if len(side) < len(other):
            canonical = side
        elif len(other) < len(side):
            canonical = other
        else:
            canonical = min(
                (side, other),
                key=lambda s: tuple(sorted(s))
            )

        unique_splits.setdefault(canonical, True)

    print("Internal unique splits:",
          len(unique_splits), flush=True)

    print("STEP 2: Calculating carrier distributions",
          flush=True)

    all_records = []
    candidate_records = []
    indexed_splits = sorted(
        unique_splits,
        key=lambda x: (len(x), tuple(sorted(x)))
    )

    for index, side_a in enumerate(indexed_splits, 1):
        split_id = f"B{index:04d}"
        side_b = all_tips - side_a

        for side_name, side in (("A", side_a), ("B", side_b)):
            assessed = side & sample_ids
            outside = sample_ids - assessed

            for hap in sorted(
                hap_ids, key=lambda x: int(x[1:])
            ):
                in_carriers = assessed & carriers[hap]
                out_carriers = outside & carriers[hap]

                n_inside = len(assessed)
                n_outside = len(outside)

                if n_inside == 0 or n_outside == 0:
                    continue

                prop_inside = len(in_carriers) / n_inside
                prop_outside = len(out_carriers) / n_outside
                coverage = (
                    len(in_carriers) / len(carriers[hap])
                )

                row = {
                    "split_id": split_id,
                    "side": side_name,
                    "haplotype_id": hap,
                    "tree_tips_on_side": len(side),
                    "assessed_on_side": n_inside,
                    "assessed_outside": n_outside,
                    "carriers_on_side": len(in_carriers),
                    "carriers_outside": len(out_carriers),
                    "carrier_proportion_on_side":
                        f"{prop_inside:.9f}",
                    "carrier_proportion_outside":
                        f"{prop_outside:.9f}",
                    "proportion_difference":
                        f"{prop_inside-prop_outside:.9f}",
                    "carrier_coverage":
                        f"{coverage:.9f}",
                    "carrier_tip_ids":
                        ";".join(sorted(in_carriers)),
                    "all_assessed_tip_ids":
                        ";".join(sorted(assessed)),
                }
                all_records.append(row)

                if (
                    n_inside >= 5
                    and n_outside >= 5
                    and len(in_carriers) >= 3
                    and prop_inside > prop_outside
                ):
                    candidate_records.append(row)

    require(all_records,
            "No branch correspondence generated")

    print("STEP 3: Selecting descriptive summaries",
          flush=True)

    best = {}
    for hap in sorted(hap_ids, key=lambda x: int(x[1:])):
        options = [
            r for r in candidate_records
            if r["haplotype_id"] == hap
        ]
        options.sort(
            key=lambda r: (
                -float(r["proportion_difference"]),
                -int(r["carriers_on_side"]),
                int(r["assessed_on_side"]),
                r["split_id"],
                r["side"]
            )
        )
        best[hap] = options[0] if options else None

    summary = []

    for h in catalogue:
        hid = h["haplotype_id"]
        b = best[hid]
        row = {
            "haplotype_id": hid,
            "amino_acid_haplotype":
                h["amino_acid_haplotype"],
            "chromosome_count": h["chromosome_count"],
            "carrier_count": len(carriers[hid]),
            "carrier_tip_ids":
                ";".join(sorted(carriers[hid])),
            "selected_split_id":
                b["split_id"] if b else "",
            "selected_side":
                b["side"] if b else "",
            "assessed_on_selected_side":
                b["assessed_on_side"] if b else "",
            "carriers_on_selected_side":
                b["carriers_on_side"] if b else "",
            "carrier_proportion_on_selected_side":
                b["carrier_proportion_on_side"] if b else "",
            "carrier_coverage_on_selected_side":
                b["carrier_coverage"] if b else "",
            "proportion_difference":
                b["proportion_difference"] if b else "",
            "selection_status":
                "descriptive_candidate" if b
                else "no_eligible_split",
        }
        summary.append(row)

    write_tsv(
        out / "H1_H26_branch_correspondence_summary.tsv",
        summary, list(summary[0])
    )
    write_tsv(
        out / "all_branch_haplotype_sides.tsv",
        all_records, list(all_records[0])
    )

    selected = sorted(
        candidate_records,
        key=lambda r: (
            int(r["haplotype_id"][1:]),
            -float(r["proportion_difference"]),
            r["split_id"],
            r["side"]
        )
    )
    write_tsv(
        out / "eligible_enriched_branch_sides.tsv",
        selected, list(all_records[0])
    )

    split_summary = []
    for i, side in enumerate(indexed_splits, 1):
        split_summary.append({
            "split_id": f"B{i:04d}",
            "side_A_total_tree_tips": len(side),
            "side_B_total_tree_tips":
                len(all_tips) - len(side),
            "side_A_assessed":
                len(side & sample_ids),
            "side_B_assessed":
                len(sample_ids - side),
        })

    write_tsv(
        out / "branch_split_inventory.tsv",
        split_summary, list(split_summary[0])
    )

    report = {
        "generated_utc":
            datetime.now(timezone.utc).isoformat(),
        "tree_tips": 1217,
        "assessed_tips": len(sample_ids),
        "unassessed_tips": 1217 - len(sample_ids),
        "protein_haplotypes": len(hap_ids),
        "unique_internal_splits": len(indexed_splits),
        "evaluated_side_haplotype_rows":
            len(all_records),
        "eligible_enriched_side_rows":
            len(candidate_records),
        "haplotypes_with_candidate_split":
            sum(b is not None for b in best.values()),
        "criteria": {
            "minimum_assessed_on_side": 5,
            "minimum_assessed_outside": 5,
            "minimum_carriers_on_side": 3,
            "require_higher_inside_proportion": True
        },
        "interpretation": [
            "Unrooted branch splits do not imply ancestry.",
            "Only 165 assessed cattle have haplotype labels.",
            "The selected split is descriptive, not a test.",
            "Selecting the maximum contrast over branches "
            "introduces selection bias.",
            "Carrier groups overlap across haplotypes.",
            "No branch-level significance is claimed."
        ]
    }

    (out / "branch_correspondence_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )

    print(json.dumps(report, indent=2), flush=True)
    print("BRANCH CORRESPONDENCE COMPLETED", flush=True)


if __name__ == "__main__":
    main()
