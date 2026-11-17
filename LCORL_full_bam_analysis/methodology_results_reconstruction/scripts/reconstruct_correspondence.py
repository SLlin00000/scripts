#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import re
import hashlib
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime, timezone

SEGMENTS = [f"S{i}" for i in range(1, 13)]

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows):
    check(bool(rows), f"No rows for {path}")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]),
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    root = args.project.resolve()
    base = root / "LCORL_full_bam_analysis"
    treebase = root / "cattle_phylogeny"
    out = args.out.resolve()

    check(out.is_dir() and not any(out.iterdir()),
          "Output directory must exist and be empty")

    files = {
        "vcf": base / "haplotype/LCORL_gene_phased_MAF01.vcf.gz",
        "evidence": base / "summary/functional_annotation/LCORL_representative_missense_MAF01_evidence.tsv",
        "segments": base / "haplotype/segment_SNP_annotation.tsv",
        "definitions": base / "haplotype/LCORL_adjacent_highLD_segments.tsv",
        "old_protein": base / "summary/functional_annotation/LCORL_joint_missense_protein_haplotypes_MAF01.tsv",
        "old_mapping": base / "summary/functional_annotation/LCORL_joint_missense_segment_mapping_MAF01.tsv",
        "metadata": treebase / "metadata/tree1217_annotation_final.tsv",
        "tree_ids": treebase / "metadata/tree1217_sample_ids.txt",
        "newick": treebase / "results/cattle1217_LCORL_NJ_tree.newick",
    }

    for name, path in files.items():
        check(path.is_file(), f"Missing input: {name}: {path}")

    print("STEP 1: Loading variant annotations", flush=True)
    evidence = read_tsv(files["evidence"])
    check(len(evidence) == 11, "Expected 11 missense variants")

    variants = []
    for r in evidence:
        m = re.fullmatch(r"([A-Z*])(\d+)([A-Z*])",
                         r["Protein_change"])
        check(m is not None, "Invalid protein change")
        variants.append({
            "pos": int(r["Genomic_position"]),
            "ref": r["Ref"],
            "alt": r["Alt"],
            "aa_pos": int(m.group(2)),
            "refaa": m.group(1),
            "altaa": m.group(3),
            "change": r["Protein_change"],
            "row": r
        })

    variants.sort(key=lambda x: x["aa_pos"])
    check(len({v["pos"] for v in variants}) == 11,
          "Duplicate missense positions")

    print("STEP 2: Checking 12 LD segments", flush=True)
    segment_rows = read_tsv(files["segments"])
    positions = defaultdict(list)
    for r in segment_rows:
        positions[r["Segment"]].append(int(r["Position"]))

    check(set(positions) == set(SEGMENTS),
          "Unexpected segment IDs")
    positions = {
        s: sorted(set(positions[s])) for s in SEGMENTS
    }
    check(sum(map(len, positions.values())) == 43,
          "Expected 43 segment markers")

    definitions = read_tsv(files["definitions"])
    check(len(definitions) == 12, "Expected 12 segments")

    for d in definitions:
        s = d["SEGMENT"]
        check(s in positions, f"Unknown segment: {s}")
        listed = [int(x) for x in d["SNP_POSITIONS"].split(",")]
        check(listed == positions[s],
              f"Segment definition mismatch: {s}")

    wanted = {v["pos"] for v in variants}
    wanted.update(p for ps in positions.values() for p in ps)

    print("STEP 3: Reading phased VCF", flush=True)
    genotypes = {}
    nucleotides = {}
    samples = None

    with gzip.open(files["vcf"], "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                check(len(samples) == 165, "Expected 165 samples")
                check(len(set(samples)) == 165,
                      "Duplicate sample IDs")
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            pos = int(fields[1])
            if pos not in wanted:
                continue

            check(fields[0] == "6", "Unexpected chromosome")
            check(pos not in genotypes,
                  f"Duplicate position: {pos}")

            nt = [fields[3], *fields[4].split(",")]
            check(len(nt) == 2 and all(len(x) == 1 for x in nt),
                  f"Non-biallelic SNP: {pos}")

            fmt = fields[8].split(":")
            check("GT" in fmt, f"No GT at {pos}")
            index = fmt.index("GT")

            calls = []
            for field in fields[9:]:
                values = field.split(":")
                check(index < len(values), f"Missing GT at {pos}")
                gt = values[index]
                check(re.fullmatch(r"[01]\|[01]", gt) is not None,
                      f"Unphased or missing GT at {pos}: {gt}")
                calls.append((int(gt[0]), int(gt[2])))

            check(len(calls) == len(samples),
                  f"Sample count mismatch at {pos}")
            genotypes[pos] = calls
            nucleotides[pos] = nt

    check(set(genotypes) == wanted,
          f"Missing VCF positions: {sorted(wanted-set(genotypes))}")

    for v in variants:
        check(nucleotides[v["pos"]] == [v["ref"], v["alt"]],
              f"REF/ALT mismatch: {v['pos']}")

    print("STEP 4: Validating tree metadata", flush=True)
    metadata = read_tsv(files["metadata"])
    check(len(metadata) == 1217, "Expected 1217 metadata rows")
    meta = {r["sample_id"]: r for r in metadata}
    check(len(meta) == 1217, "Duplicate metadata IDs")

    with open(files["tree_ids"]) as f:
        tree_ids = [
            line.strip().split()[0]
            for line in f if line.strip()
        ]

    check(len(tree_ids) == 1217 and len(set(tree_ids)) == 1217,
          "Invalid tree ID list")
    check(set(tree_ids) == set(meta),
          "Tree list and metadata mismatch")
    check(set(samples) <= set(tree_ids),
          "Phased samples absent from tree")

    from Bio import Phylo
    tree = Phylo.read(str(files["newick"]), "newick")
    tips = [t.name for t in tree.get_terminals()]
    check(len(tips) == 1217 and len(set(tips)) == 1217,
          "Invalid Newick tips")
    check(set(tips) == set(tree_ids),
          "Newick and tree IDs mismatch")

    print("STEP 5: Reconstructing chromosome haplotypes",
          flush=True)

    chromosomes = []
    protein_counts = Counter()
    joint_counts = Counter()
    sample_patterns = defaultdict(list)
    descriptions = {}

    for i, sample in enumerate(samples):
        for copy in (0, 1):
            aa = []
            changes = []
            allele_codes = []
            nucleotide_codes = []

            for v in variants:
                code = genotypes[v["pos"]][i][copy]
                aa.append(
                    f"{v['aa_pos']}"
                    f"{v['altaa'] if code else v['refaa']}"
                )
                allele_codes.append(str(code))
                nucleotide_codes.append(
                    v["alt"] if code else v["ref"]
                )
                if code:
                    changes.append(v["change"])

            pattern = "-".join(aa)
            changes_text = (
                ",".join(changes) if changes
                else "Reference_at_all_11"
            )

            segment_patterns = {}
            for s in SEGMENTS:
                segment_patterns[s] = "-".join(
                    nucleotides[p][genotypes[p][i][copy]]
                    for p in positions[s]
                )

            key = (pattern, changes_text)
            protein_counts[key] += 1
            sample_patterns[sample].append(pattern)

            joint_counts[(
                pattern, changes_text,
                segment_patterns["S1"],
                segment_patterns["S5"],
                segment_patterns["S11"]
            )] += 1

            descriptions[pattern] = (
                "|".join(allele_codes),
                "|".join(nucleotide_codes)
            )

            row = {
                "sample_id": sample,
                "chromosome_copy": copy,
                "amino_acid_haplotype": pattern,
                "ALT_missense_changes": changes_text,
                "allele_codes_in_protein_order":
                    "|".join(allele_codes),
                "nucleotide_alleles_in_protein_order":
                    "|".join(nucleotide_codes),
            }
            row.update({
                f"{s}_haplotype": segment_patterns[s]
                for s in SEGMENTS
            })
            row.update(meta[sample])
            chromosomes.append(row)

    check(len(chromosomes) == 330,
          "Expected 330 chromosome copies")
    check(all(len(x) == 2 for x in sample_patterns.values()),
          "Incomplete individual haplotypes")

    ranking = sorted(
        protein_counts.items(),
        key=lambda x: (-x[1], x[0])
    )
    check(len(ranking) == 26, "Expected 26 haplotypes")
    check(sum(protein_counts.values()) == 330,
          "Haplotype count does not sum to 330")

    print("STEP 6: Comparing original results", flush=True)
    original = read_tsv(files["old_protein"])
    check(len(original) == 26, "Original catalogue invalid")

    hap_ids = {}
    for rank, (key, count) in enumerate(ranking, 1):
        pattern, changes = key
        hap_ids[pattern] = f"H{rank}"
        old = original[rank-1]

        check(
            int(old["Rank"]) == rank
            and old["Amino_acid_haplotype"] == pattern
            and old["ALT_missense_changes"] == changes
            and int(old["Count"]) == count
            and abs(float(old["Frequency"])-count/330) <= 5.1e-7,
            f"Original protein result mismatch: H{rank}"
        )

    old_joint = {}
    for r in read_tsv(files["old_mapping"]):
        key = (
            r["Amino_acid_haplotype"],
            r["ALT_missense_changes"],
            r["S1_haplotype"],
            r["S5_haplotype"],
            r["S11_haplotype"]
        )
        check(key not in old_joint,
              "Duplicate original mapping row")
        old_joint[key] = int(r["Count"])
        check(
            abs(float(r["Frequency"])-int(r["Count"])/330)
            <= 5.1e-7,
            "Original mapping frequency mismatch"
        )

    check(old_joint == dict(joint_counts),
          "S1/S5/S11 joint mapping differs from original")

    print("STEP 7: Building output tables", flush=True)
    for r in chromosomes:
        r["protein_haplotype_id"] = hap_ids[
            r["amino_acid_haplotype"]
        ]

    catalogue = []
    for rank, (key, count) in enumerate(ranking, 1):
        pattern, changes = key
        hid = f"H{rank}"
        carriers = {
            r["sample_id"] for r in chromosomes
            if r["protein_haplotype_id"] == hid
        }
        row = {
            "haplotype_id": hid,
            "rank_original": rank,
            "amino_acid_haplotype": pattern,
            "ALT_missense_changes": changes,
            "allele_codes_in_protein_order":
                descriptions[pattern][0],
            "nucleotide_alleles_in_protein_order":
                descriptions[pattern][1],
            "chromosome_count": count,
            "frequency": f"{count/330:.9f}",
            "carrier_individuals": len(carriers),
        }
        for v, state in zip(variants, pattern.split("-")):
            row[f"AA_{v['aa_pos']}"] = state
        catalogue.append(row)

    diplotypes = {}
    for sample in samples:
        first, second = [
            hap_ids[p] for p in sample_patterns[sample]
        ]
        ordered = sorted(
            (first, second),
            key=lambda x: int(x[1:])
        )
        diplotypes[sample] = {
            "haplotype_copy0": first,
            "haplotype_copy1": second,
            "diplotype": "/".join(ordered),
            "haplotype_homozygous": first == second,
        }

    individuals = [
        {"sample_id": s, **diplotypes[s], **meta[s]}
        for s in sorted(samples)
    ]

    tree_rows = []
    for sample in tips:
        if sample in diplotypes:
            row = {
                **meta[sample],
                "haplotype_assessment": "assessed",
                **diplotypes[sample]
            }
        else:
            row = {
                **meta[sample],
                "haplotype_assessment": "not_assessed",
                "haplotype_copy0": "",
                "haplotype_copy1": "",
                "diplotype": "",
                "haplotype_homozygous": ""
            }
        tree_rows.append(row)

    all_segments = []
    for s in SEGMENTS:
        c = Counter(
            r[f"{s}_haplotype"] for r in chromosomes
        )
        for pattern, count in sorted(
            c.items(), key=lambda x: (-x[1], x[0])
        ):
            all_segments.append({
                "segment_id": s,
                "segment_haplotype": pattern,
                "chromosome_count": count,
                "frequency": f"{count/330:.9f}",
                "marker_positions":
                    ",".join(map(str, positions[s]))
            })

    joint_rows = []
    for key, count in sorted(
        joint_counts.items(),
        key=lambda x: (-x[1], x[0])
    ):
        joint_rows.append({
            "haplotype_id": hap_ids[key[0]],
            "amino_acid_haplotype": key[0],
            "ALT_missense_changes": key[1],
            "S1_haplotype": key[2],
            "S5_haplotype": key[3],
            "S11_haplotype": key[4],
            "Count": count,
            "Frequency": f"{count/330:.9f}",
        })

    print("STEP 8: Saving validated results", flush=True)
    write_tsv(
        out / "variant_catalogue.tsv",
        [{"protein_order": i, **v["row"]}
         for i, v in enumerate(variants, 1)]
    )
    write_tsv(out / "protein_haplotype_catalogue.tsv",
              catalogue)
    write_tsv(out / "chromosome_haplotype_assignments.tsv",
              chromosomes)
    write_tsv(out / "cattle_diplotype_assignments.tsv",
              individuals)
    write_tsv(out / "tree1217_haplotype_annotations.tsv",
              tree_rows)
    write_tsv(out / "all12_segment_haplotype_frequencies.tsv",
              all_segments)
    write_tsv(out / "protein_S1_S5_S11_joint_mapping.tsv",
              joint_rows)
    write_tsv(out / "LD_segment_definitions.tsv",
              definitions)

    write_tsv(
        out / "input_file_checksums.tsv",
        [{"input_name": name, "path": str(path),
          "sha256": sha256(path)}
         for name, path in files.items()]
    )

    report = {
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "phased_samples": len(samples),
        "chromosome_copies": len(chromosomes),
        "protein_haplotypes": len(ranking),
        "missense_variants": len(variants),
        "LD_segments": len(positions),
        "LD_segment_markers":
            sum(map(len, positions.values())),
        "newick_tips": len(tips),
        "phased_cattle_in_tree":
            len(set(samples) & set(tips)),
        "tree_tips_not_assessed":
            len(tips) - len(samples),
        "original_protein_counts_match": True,
        "original_S1_S5_S11_mapping_match": True,
        "independent_phasing_accuracy_validation": False,
        "notes": [
            "H1-H26 IDs preserve original frequency ranks.",
            "Copy 0/1 is not maternal/paternal assignment.",
            "Non-phased tree tips are not_assessed.",
            "S1-S12 mapping extends the original S1/S5/S11 analysis.",
            "Breed and ancestry labels are not inferred."
        ]
    }

    (out / "validation_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )

    print(json.dumps(report, indent=2), flush=True)
    print("VALIDATION PASSED", flush=True)

if __name__ == "__main__":
    main()
