#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent

COMMON_VCF = ROOT / "haplotype/LCORL_gene_phased_MAF01.vcf.gz"
CDS_VCF = ROOT / "haplotype/LCORL_common_CDS_MAF01_alltranscripts.vcf.gz"

MISSENSE_TSV = (
    ROOT
    / "summary/functional_annotation/"
      "LCORL_representative_missense_MAF01_evidence.tsv"
)

HAPLOTYPE_TSV = (
    ROOT
    / "summary/functional_annotation/"
      "LCORL_joint_missense_protein_haplotypes_MAF01.tsv"
)

OUTPUT_DIR = ROOT / "figures"
OUTPUT_PNG = OUTPUT_DIR / "LCORL_Q9774_summary_from_data_MAF01.png"
OUTPUT_PDF = OUTPUT_DIR / "LCORL_Q9774_summary_from_data_MAF01.pdf"
OUTPUT_STATS = OUTPUT_DIR / "LCORL_Q9774_summary_plot_statistics_MAF01.tsv"

PROTEIN_LENGTH = 1904
DUF4553_START = 1442
DUF4553_END = 1897

# Previously verified by the exon/intron annotation step.
EXON_NON_CDS_COUNT = 5


def require_file(path):
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_vcf_positions(path):
    """Return unique (chromosome, position) records from a VCF."""
    positions = set()

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            chrom = fields[0]
            position = int(fields[1])
            positions.add((chrom, position))

    return positions


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_columns(rows, required, filename):
    if not rows:
        raise ValueError(f"No rows found in {filename}")

    missing = [column for column in required if column not in rows[0]]

    if missing:
        available = ", ".join(rows[0].keys())
        raise ValueError(
            f"{filename} is missing columns: {missing}\n"
            f"Available columns: {available}"
        )


def protein_position(change):
    """Convert K1503Q to 1503."""
    match = re.search(r"(\d+)", change)

    if not match:
        raise ValueError(f"Cannot obtain protein position from: {change}")

    return int(match.group(1))


def is_inside_duf(row):
    text = row.get("DUF4553_status", "").lower()
    return "within" in text or "inside" in text


def prepare_data():
    for path in [COMMON_VCF, CDS_VCF, MISSENSE_TSV, HAPLOTYPE_TSV]:
        require_file(path)

    common_positions = read_vcf_positions(COMMON_VCF)
    cds_positions = read_vcf_positions(CDS_VCF)

    total_common = len(common_positions)
    cds_count = len(cds_positions)
    intronic_count = total_common - cds_count - EXON_NON_CDS_COUNT

    if intronic_count < 0:
        raise ValueError(
            "Calculated intronic count is negative. Check input VCF files."
        )

    missense = read_tsv(MISSENSE_TSV)
    require_columns(
        missense,
        ["Protein_change", "MAF", "SIFT_score"],
        MISSENSE_TSV.name,
    )

    for row in missense:
        row["Protein_position_numeric"] = protein_position(
            row["Protein_change"]
        )
        row["MAF_numeric"] = float(row["MAF"])
        row["SIFT_numeric"] = float(row["SIFT_score"])

    missense.sort(key=lambda row: row["Protein_position_numeric"])

    haplotypes = read_tsv(HAPLOTYPE_TSV)
    require_columns(
        haplotypes,
        [
            "Rank",
            "Amino_acid_haplotype",
            "ALT_missense_changes",
            "Count",
            "Frequency",
        ],
        HAPLOTYPE_TSV.name,
    )

    haplotypes.sort(key=lambda row: int(row["Rank"]))

    haplotype_total = sum(int(row["Count"]) for row in haplotypes)

    annotation_counts = {
        "Intronic": intronic_count,
        "CDS": cds_count,
        "Exon non-CDS": EXON_NON_CDS_COUNT,
    }

    return (
        annotation_counts,
        total_common,
        missense,
        haplotypes,
        haplotype_total,
    )


def write_statistics(
    annotation_counts,
    total_common,
    missense,
    haplotypes,
    haplotype_total,
):
    with open(OUTPUT_STATS, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Section", "Category", "Count_or_value"])

        for category, count in annotation_counts.items():
            writer.writerow(["Variant_annotation", category, count])

        writer.writerow(
            ["Variant_annotation", "Total_common_variants", total_common]
        )

        for row in missense:
            writer.writerow(
                [
                    "Missense_variant",
                    row["Protein_change"],
                    f"MAF={row['MAF_numeric']:.6f};"
                    f"SIFT={row['SIFT_numeric']:.6f}",
                ]
            )

        for row in haplotypes:
            writer.writerow(
                [
                    "Protein_haplotype",
                    row["Amino_acid_haplotype"],
                    f"count={row['Count']};frequency={row['Frequency']}",
                ]
            )

        writer.writerow(
            ["Protein_haplotype", "Total_chromosome_copies", haplotype_total]
        )


def make_figure(
    annotation_counts,
    total_common,
    missense,
    haplotypes,
    haplotype_total,
):
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(18, 11.5))
    fig.subplots_adjust(
        left=0.09,
        right=0.97,
        bottom=0.14,
        top=0.88,
        wspace=0.55,
        hspace=0.52,
    )

    # Panel A: exon/intron/CDS distribution
    ax = axes[0, 0]

    classes = list(annotation_counts.keys())
    counts = np.array(list(annotation_counts.values()))
    colours = ["#4C78A8", "#E45756", "#72B7B2"]
    y = np.arange(len(classes))

    ax.barh(y, counts, color=colours, height=0.62)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(classes)
    ax.invert_yaxis()
    ax.set_xlabel("Number of variants (log scale)")
    ax.set_title(
        f"A  Functional location of variants with MAF ≥ 0.01 (n={total_common})",
        loc="left",
    )
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_xlim(0.8, max(counts) * 2)

    for yi, value in zip(y, counts):
        percentage = 100 * value / total_common
        ax.text(
            value * 1.10,
            yi,
            f"{value} ({percentage:.1f}%)",
            va="center",
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B: representative LCORL protein
    ax = axes[0, 1]

    ax.set_xlim(0, PROTEIN_LENGTH + 40)
    ax.set_ylim(-0.18, 1.48)

    ax.hlines(
        0.48,
        1,
        PROTEIN_LENGTH,
        color="#59636E",
        linewidth=5,
        zorder=1,
    )

    ax.add_patch(
        plt.Rectangle(
            (DUF4553_START, 0.34),
            DUF4553_END - DUF4553_START,
            0.28,
            facecolor="#8B6BB7",
            edgecolor="#5B3E85",
            linewidth=1.2,
            zorder=2,
        )
    )

    ax.text(
        (DUF4553_START + DUF4553_END) / 2,
        0.48,
        "DUF4553\n1442–1897 aa",
        ha="center",
        va="center",
        color="white",
        fontsize=8.5,
        fontweight="bold",
        zorder=3,
    )

    label_levels = [0.90, 1.20]

    for index, row in enumerate(missense):
        position = row["Protein_position_numeric"]
        change = row["Protein_change"]
        inside_domain = is_inside_duf(row)

        colour = "#C23B3B" if inside_domain else "#2F6690"
        label_y = label_levels[index % 2]

        ax.vlines(
            position,
            0.49,
            label_y - 0.06,
            color=colour,
            linewidth=1.4,
        )
        ax.scatter(
            position,
            0.49,
            color=colour,
            edgecolor="white",
            linewidth=0.7,
            s=52 if inside_domain else 38,
            zorder=4,
        )
        ax.text(
            position,
            label_y,
            change,
            ha="center",
            va="bottom",
            color=colour,
            fontsize=8.2,
            rotation=55,
            rotation_mode="anchor",
            fontweight="bold" if inside_domain else "normal",
        )

    ax.set_xticks([0, 500, 1000, 1500, PROTEIN_LENGTH])
    ax.set_xlabel("Protein position (aa)")
    ax.set_yticks([])
    ax.set_title(
        "B  Representative LCORL protein and 11 missense variants",
        loc="left",
    )

    ax.text(
        1,
        0.20,
        "1 aa",
        ha="left",
        va="top",
        color="#59636E",
    )
    ax.text(
        PROTEIN_LENGTH,
        0.20,
        f"{PROTEIN_LENGTH} aa",
        ha="right",
        va="top",
        color="#59636E",
    )

    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Panel C: missense MAF and SIFT
    ax = axes[1, 0]

    names = [row["Protein_change"] for row in missense]
    maf = np.array([row["MAF_numeric"] for row in missense])
    sift = np.array([row["SIFT_numeric"] for row in missense])
    inside = [is_inside_duf(row) for row in missense]

    x = np.arange(len(names))
    colours = [
        "#C23B3B" if domain_status else "#4C78A8"
        for domain_status in inside
    ]

    bars = ax.bar(x, maf, color=colours, width=0.68)

    ax.axhline(
        0.01,
        color="#777777",
        linestyle="--",
        linewidth=1,
        label="MAF threshold = 0.01",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=42, ha="right", fontsize=8.5)
    ax.set_ylabel("Minor allele frequency")
    ax.set_ylim(0, max(0.62, max(maf) + 0.13))
    ax.set_title(
        "C  Missense allele frequency and SIFT evidence",
        loc="left",
    )
    ax.grid(axis="y", color="#E1E1E1", linewidth=0.7)
    ax.set_axisbelow(True)

    for rect, maf_value, sift_value in zip(bars, maf, sift):
        ax.text(
            rect.get_x() + rect.get_width() / 2,
            maf_value + 0.012,
            f"MAF {maf_value:.3f}\nSIFT {sift_value:g}",
            ha="center",
            va="bottom",
            fontsize=7.2,
            linespacing=1.15,
        )

    ax.legend(frameon=False, loc="upper right", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel D: major protein haplotypes
    ax = axes[1, 1]

    top_three = haplotypes[:3]
    top_counts = [int(row["Count"]) for row in top_three]
    other_count = haplotype_total - sum(top_counts)

    hap_labels = []
    for row in top_three:
        changes = row["ALT_missense_changes"]

        if changes == "I562V,A762P,S913G,E1019K,K1503Q":
            hap_labels.append(
                "5-variant haplotype\n"
                "I562V, A762P, S913G, E1019K, K1503Q"
            )
        elif changes == "A762P,S913G,K1503Q":
            hap_labels.append(
                "3-variant haplotype\n"
                "A762P, S913G, K1503Q"
            )
        elif changes == "Reference_at_all_11":
            hap_labels.append(
                "Reference haplotype\n"
                "reference amino acid at all 11 sites"
            )
        else:
            hap_labels.append(changes.replace(",", ", "))

    hap_labels.append(
        f"Other {max(0, len(haplotypes) - 3)} haplotypes"
    )
    plot_counts = np.array(top_counts + [other_count])

    colours = ["#5B8E7D", "#E9A03B", "#8D99AE", "#D4D4D4"]
    y = np.arange(len(hap_labels))

    ax.barh(y, plot_counts, color=colours, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels(hap_labels, fontsize=8.5)
    ax.tick_params(axis="y", pad=9)

    for tick_label in ax.get_yticklabels():
        tick_label.set_linespacing(1.20)

    ax.invert_yaxis()
    ax.set_xlabel(
        f"Phased chromosome copies (total={haplotype_total})"
    )
    ax.set_title("D  Major phased protein haplotypes", loc="left")
    ax.grid(axis="x", color="#E1E1E1", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max(plot_counts) * 1.38)

    for yi, value in zip(y, plot_counts):
        percentage = 100 * value / haplotype_total
        ax.text(
            value + 2,
            yi,
            f"{value}/{haplotype_total} ({percentage:.1f}%)",
            va="center",
            fontsize=9,
        )


    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle(
        "LCORL functional annotation and phased protein haplotypes",
        fontsize=18,
        fontweight="bold",
        y=0.97,
    )

    fig.text(
        0.5,
        0.925,
        "Q9774 discovery cohort: 165 cattle / "
        "330 phased chromosome copies; ARS-UCD1.2",
        ha="center",
        fontsize=11,
        color="#444444",
    )

    return fig


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    (
        annotation_counts,
        total_common,
        missense,
        haplotypes,
        haplotype_total,
    ) = prepare_data()

    if total_common != 1011:
        print(
            f"WARNING: expected 1011 variants with MAF ≥ 0.01, "
            f"but found {total_common}",
            file=sys.stderr,
        )

    if len(missense) != 11:
        print(
            f"WARNING: expected 11 representative missense variants, "
            f"but found {len(missense)}",
            file=sys.stderr,
        )

    if haplotype_total != 330:
        print(
            f"WARNING: expected 330 phased chromosome copies, "
            f"but found {haplotype_total}",
            file=sys.stderr,
        )

    write_statistics(
        annotation_counts,
        total_common,
        missense,
        haplotypes,
        haplotype_total,
    )

    fig = make_figure(
        annotation_counts,
        total_common,
        missense,
        haplotypes,
        haplotype_total,
    )

    fig.savefig(
        OUTPUT_PNG,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        OUTPUT_PDF,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    print("Input data checked successfully")
    print(f"Common variants: {total_common}")
    print(f"Annotation counts: {annotation_counts}")
    print(f"Missense variants: {len(missense)}")
    print(f"Protein haplotypes: {len(haplotypes)}")
    print(f"Chromosome copies: {haplotype_total}")
    print(f"Created: {OUTPUT_PNG}")
    print(f"Created: {OUTPUT_PDF}")
    print(f"Created: {OUTPUT_STATS}")


if __name__ == "__main__":
    main()

