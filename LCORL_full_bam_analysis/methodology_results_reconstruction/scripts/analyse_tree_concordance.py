#!/usr/bin/env python3
import argparse
import csv
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime, timezone

from Bio import Phylo


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def mean(values):
    return sum(values) / len(values) if values else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=7621)
    args = parser.parse_args()

    require(args.permutations > 0, "Invalid permutation count")

    project = args.project.resolve()
    base = project / "LCORL_full_bam_analysis"
    original = (
        base / "methodology_results_reconstruction"
        / "results/run_30012993"
    )
    treefile = (
        project / "cattle_phylogeny/results/"
        "cattle1217_LCORL_NJ_tree.newick"
    )
    out = args.out.resolve()

    require(out.is_dir() and not any(out.iterdir()),
            "Output directory must be empty")

    require(
        (original / "COMPLETED_VALIDATED.txt").read_text().strip()
        == "PASS",
        "Original reconstruction not validated"
    )

    print("STEP 1: Loading original mapping", flush=True)

    individuals = load_tsv(
        original / "cattle_diplotype_assignments.tsv"
    )
    haplotypes = load_tsv(
        original / "protein_haplotype_catalogue.tsv"
    )

    require(len(individuals) == 165, "Expected 165 cattle")
    require(len(haplotypes) == 26, "Expected 26 haplotypes")

    sample_map = {r["sample_id"]: r for r in individuals}
    require(len(sample_map) == 165, "Duplicate sample IDs")

    hap_map = {
        r["haplotype_id"]: r for r in haplotypes
    }
    require(len(hap_map) == 26, "Duplicate haplotype IDs")

    print("STEP 2: Reading original NJ tree", flush=True)

    tree = Phylo.read(str(treefile), "newick")
    terminals = tree.get_terminals()
    tips = [x.name for x in terminals]

    require(len(tips) == 1217, "Expected 1217 tips")
    require(len(set(tips)) == 1217, "Duplicate tree tips")
    require(set(sample_map) <= set(tips),
            "Phased cattle missing from tree")

    branches = [
        c.branch_length
        for c in tree.find_clades()
        if c is not tree.root
    ]

    require(
        all(
            x is not None and math.isfinite(x)
            for x in branches
        ),
        "Missing or nonfinite branch lengths"
    )

    negative_branches = sum(x < 0 for x in branches)

    print(
        "Negative branch lengths:",
        negative_branches,
        flush=True
    )

    ids = sorted(sample_map)
    labels = [sample_map[s]["diplotype"] for s in ids]

    counts = Counter(labels)

    require(len(counts) == 47,
            "Expected 47 diplotypes")

    print("STEP 3: Creating H1-H26 table", flush=True)

    correspondence = []
    for hid in sorted(hap_map, key=lambda s: int(s[1:])):
        h = hap_map[hid]
        carriers = [
            s for s in ids
            if hid in (
                sample_map[s]["haplotype_copy0"],
                sample_map[s]["haplotype_copy1"]
            )
        ]

        diplotype_counts = Counter(
            sample_map[s]["diplotype"] for s in carriers
        )

        require(
            len(carriers) == int(h["carrier_individuals"]),
            f"Carrier count mismatch for {hid}"
        )

        correspondence.append({
            "haplotype_id": hid,
            "amino_acid_haplotype":
                h["amino_acid_haplotype"],
            "ALT_missense_changes":
                h["ALT_missense_changes"],
            "chromosome_count":
                h["chromosome_count"],
            "chromosome_frequency":
                h["frequency"],
            "carrier_count": len(carriers),
            "carrier_frequency":
                f"{len(carriers)/165:.9f}",
            "diplotype_combinations": ";".join(
                f"{k}:{v}" for k, v in
                sorted(diplotype_counts.items())
            ),
            "tree_tip_ids": ";".join(carriers),
        })

    save_tsv(
        out / "H1_H26_phylogenetic_correspondence.tsv",
        correspondence,
        list(correspondence[0])
    )

    mapping = [
        {
            "tree_tip_id": s,
            "sample_id": s,
            "haplotype_copy0":
                sample_map[s]["haplotype_copy0"],
            "haplotype_copy1":
                sample_map[s]["haplotype_copy1"],
            "diplotype":
                sample_map[s]["diplotype"],
        }
        for s in ids
    ]

    save_tsv(
        out / "individual_tree_mapping.tsv",
        mapping,
        list(mapping[0])
    )

    print("STEP 4: Calculating patristic distances",
          flush=True)

    tip_lookup = {tip.name: tip for tip in terminals}
    pairs = []
    pair_rows = []
    within_by_group = defaultdict(list)

    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = ids[i], ids[j]
            d = tree.distance(
                tip_lookup[a], tip_lookup[b]
            )

            require(math.isfinite(d),
                    "Nonfinite patristic distance")

            same = labels[i] == labels[j]

            pairs.append((i, j, d))

            if same:
                within_by_group[labels[i]].append(d)

            pair_rows.append({
                "sample_1": a,
                "sample_2": b,
                "diplotype_1": labels[i],
                "diplotype_2": labels[j],
                "same_diplotype": int(same),
                "patristic_distance": f"{d:.12g}",
            })

    require(len(pairs) == 13530,
            "Incorrect pair count")

    save_tsv(
        out / "pairwise_patristic_distances.tsv",
        pair_rows, list(pair_rows[0])
    )

    group_summary = []

    for group, n in sorted(
        counts.items(),
        key=lambda x: (-x[1], x[0])
    ):
        distances = within_by_group.get(group, [])
        expected_pairs = n * (n - 1) // 2

        require(len(distances) == expected_pairs,
                f"Pair count mismatch: {group}")

        group_summary.append({
            "diplotype": group,
            "n_cattle": n,
            "within_pairs": len(distances),
            "mean_within_distance":
                f"{mean(distances):.12g}"
                if distances else "",
            "min_within_distance":
                f"{min(distances):.12g}"
                if distances else "",
            "max_within_distance":
                f"{max(distances):.12g}"
                if distances else "",
        })

    save_tsv(
        out / "diplotype_distance_summary.tsv",
        group_summary, list(group_summary[0])
    )

    print("STEP 5: Computing observed statistic",
          flush=True)

    total_sum = sum(p[2] for p in pairs)
    total_pairs = len(pairs)

    def delta_for(label_vector):
        within_sum = 0.0
        within_n = 0

        for i, j, distance in pairs:
            if label_vector[i] == label_vector[j]:
                within_sum += distance
                within_n += 1

        require(
            0 < within_n < total_pairs,
            "No valid within/between comparison"
        )

        between_sum = total_sum - within_sum
        between_n = total_pairs - within_n

        within_mean = within_sum / within_n
        between_mean = between_sum / between_n

        return (
            within_mean - between_mean,
            within_mean,
            between_mean,
            within_n,
            between_n
        )

    observed = delta_for(labels)

    print(
        "Observed delta:", observed[0],
        flush=True
    )

    print("STEP 6: Permutation test", flush=True)

    rng = random.Random(args.seed)
    permuted = labels.copy()
    null_values = []
    extreme = 0

    for k in range(1, args.permutations + 1):
        rng.shuffle(permuted)
        statistic = delta_for(permuted)[0]
        null_values.append(statistic)

        if statistic <= observed[0]:
            extreme += 1

    pvalue = (
        (extreme + 1) /
        (args.permutations + 1)
    )

    null_rows = [
        {
            "permutation": i,
            "delta_within_minus_between": f"{x:.12g}"
        }
        for i, x in enumerate(null_values, 1)
    ]

    save_tsv(
        out / "permutation_null_distribution.tsv",
        null_rows, list(null_rows[0])
    )

    print("STEP 7: Writing report", flush=True)

    report = {
        "generated_utc":
            datetime.now(timezone.utc).isoformat(),
        "original_tree_tips": len(tips),
        "mapped_cattle": len(ids),
        "unassessed_tree_tips": len(tips) - len(ids),
        "protein_haplotypes": len(haplotypes),
        "distinct_diplotypes": len(counts),
        "pairwise_distances": total_pairs,
        "within_pairs": observed[3],
        "between_pairs": observed[4],
        "mean_within_distance": observed[1],
        "mean_between_distance": observed[2],
        "observed_delta": observed[0],
        "permutations": args.permutations,
        "random_seed": args.seed,
        "one_sided_pvalue": pvalue,
        "null_mean": mean(null_values),
        "negative_tree_branches": negative_branches,
        "null_hypothesis":
            "Diplotype labels are exchangeable among "
            "the 165 fixed tree tips.",
        "alternative":
            "Within-diplotype tree distances are "
            "smaller than between-diplotype distances.",
        "limitations": [
            "The tree and haplotypes derive from the "
            "same LCORL genomic region.",
            "This is an exploratory concordance test, "
            "not an independent validation.",
            "Singleton diplotypes have no within-group pairs.",
            "No breed or population structure adjustment.",
            "The tree is fixed; topology and branch-length "
            "uncertainty are not propagated.",
            "Pairwise observations are dependent; "
            "permutation is performed at individual-label level."
        ]
    }

    (out / "tree_distance_concordance_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )

    print(json.dumps(report, indent=2), flush=True)
    print("CONCORDANCE ANALYSIS COMPLETED", flush=True)


if __name__ == "__main__":
    main()
