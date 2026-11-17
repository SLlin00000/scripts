#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, FuncNorm
import numpy as np


ROOT = Path(__file__).resolve().parent.parent

VCF = ROOT / "haplotype/LCORL_gene_phased_MAF01.vcf.gz"
EVIDENCE = (
    ROOT / "summary/functional_annotation/"
    "LCORL_representative_missense_MAF01_evidence.tsv"
)
OLD_MATRIX = ROOT / "figures/LCORL_six_missense_LD_r2.tsv"

OUTPUT_MATRIX = ROOT / "figures/LCORL_11_missense_LD_r2_MAF01.tsv"
OUTPUT_PNG = ROOT / "figures/LCORL_11_missense_LD_heatmap_MAF01.png"
OUTPUT_PDF = ROOT / "figures/LCORL_11_missense_LD_heatmap_MAF01.pdf"


def protein_position(change):
    match = re.fullmatch(r"[A-Z*](\d+)[A-Z*]", change)
    if not match:
        raise ValueError(f"Cannot parse protein change: {change}")
    return int(match.group(1))


def read_variants():
    with EVIDENCE.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if len(rows) != 11:
        raise ValueError(
            f"Expected 11 missense variants, found {len(rows)}"
        )

    variants = []
    seen_positions = set()

    for row in rows:
        position = int(row["Genomic_position"])

        if position in seen_positions:
            raise ValueError(f"Duplicate position: {position}")
        seen_positions.add(position)

        variants.append({
            "position": position,
            "ref": row["Ref"],
            "alt": row["Alt"],
            "change": row["Protein_change"],
            "protein_position": protein_position(row["Protein_change"]),
        })

    variants.sort(key=lambda x: x["protein_position"])
    return variants


def read_phased_alt_alleles(variants):
    wanted = {variant["position"]: variant for variant in variants}
    haplotypes = {}
    sample_count = None

    with gzip.open(VCF, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                sample_count = len(fields[9:])
                continue

            fields = line.rstrip().split("\t")
            position = int(fields[1])

            if position not in wanted:
                continue

            variant = wanted[position]
            ref = fields[3]
            alt = fields[4]

            if ref != variant["ref"] or alt != variant["alt"]:
                raise ValueError(
                    f"Allele mismatch at {position}: "
                    f"VCF={ref}>{alt}, evidence="
                    f"{variant['ref']}>{variant['alt']}"
                )

            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise ValueError(f"GT missing at {position}")

            gt_index = format_fields.index("GT")
            allele_vector = []

            for sample_field in fields[9:]:
                values = sample_field.split(":")
                genotype = values[gt_index]

                if "|" not in genotype:
                    raise ValueError(
                        f"Unphased or missing genotype at {position}: "
                        f"{genotype}"
                    )

                codes = genotype.split("|")
                if len(codes) != 2 or any(code not in {"0", "1"} for code in codes):
                    raise ValueError(
                        f"Unsupported genotype at {position}: {genotype}"
                    )

                allele_vector.extend(int(code) for code in codes)

            haplotypes[position] = np.asarray(
                allele_vector,
                dtype=float,
            )

    if sample_count is None:
        raise ValueError("VCF sample header was not found")

    missing = sorted(set(wanted) - set(haplotypes))
    if missing:
        raise ValueError(f"Variants missing from VCF: {missing}")

    lengths = {len(vector) for vector in haplotypes.values()}
    if lengths != {2 * sample_count}:
        raise ValueError(
            f"Inconsistent chromosome-copy counts: {sorted(lengths)}"
        )

    return haplotypes, sample_count


def pairwise_r2(x, y):
    p_a = float(np.mean(x))
    p_b = float(np.mean(y))
    p_ab = float(np.mean(x * y))

    denominator = p_a * (1.0 - p_a) * p_b * (1.0 - p_b)

    if denominator <= 0:
        return np.nan

    d_value = p_ab - p_a * p_b
    return (d_value * d_value) / denominator


def calculate_matrix(variants, haplotypes):
    number = len(variants)
    matrix = np.empty((number, number), dtype=float)

    for i, first in enumerate(variants):
        for j, second in enumerate(variants):
            if i == j:
                matrix[i, j] = 1.0
            else:
                matrix[i, j] = pairwise_r2(
                    haplotypes[first["position"]],
                    haplotypes[second["position"]],
                )

    return matrix


def validate_against_old_matrix(variants, matrix):
    if not OLD_MATRIX.is_file():
        raise FileNotFoundError(
            f"Old validation matrix missing: {OLD_MATRIX}"
        )

    with OLD_MATRIX.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    old_values = {}

    for row in rows:
        row_name = row["Variant"]
        for column_name, value in row.items():
            if column_name == "Variant":
                continue
            old_values[(row_name, column_name)] = float(value)

    index = {
        variant["change"]: i
        for i, variant in enumerate(variants)
    }

    differences = []

    for (first, second), old_value in old_values.items():
        if first not in index or second not in index:
            continue

        new_value = matrix[index[first], index[second]]
        differences.append(abs(new_value - old_value))

    if not differences:
        raise ValueError(
            "No overlapping entries found for old-matrix validation"
        )

    maximum_difference = max(differences)

    print(
        "Maximum absolute difference from old six-site matrix:",
        f"{maximum_difference:.10f}",
    )

    if maximum_difference > 2e-6:
        raise RuntimeError(
            "New calculation does not reproduce the old six-site matrix; "
            "heatmap was not generated."
        )


def write_matrix(variants, matrix):
    labels = [variant["change"] for variant in variants]

    with OUTPUT_MATRIX.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Variant"] + labels)

        for label, row in zip(labels, matrix):
            writer.writerow(
                [label] + [f"{value:.6f}" for value in row]
            )


def plot_heatmap(variants, matrix):
    labels = [variant["change"] for variant in variants]

    # Piecewise nonlinear scale expands closely spaced low-r2 values.
    # Exact r2 values remain printed inside each cell.
    data_anchors = np.array([
        0.000, 0.001, 0.002, 0.003,
        0.006, 0.007, 0.010, 0.013, 0.014,
        0.018, 0.025, 0.050, 0.100, 0.200,
        0.500, 0.750, 0.900, 1.000,
    ])
    colour_anchors = np.linspace(
        0.0, 1.0, len(data_anchors)
    )

    def colour_forward(values):
        return np.interp(
            values, data_anchors, colour_anchors
        )

    def colour_inverse(values):
        return np.interp(
            values, colour_anchors, data_anchors
        )

    colour_norm = FuncNorm(
        (colour_forward, colour_inverse),
        vmin=0.0,
        vmax=1.0,
    )

    colour_map = LinearSegmentedColormap.from_list(
        "expanded_low_r2_red_to_blue",
        [
            (0.00, "#67001F"),
            (0.10, "#B2182B"),
            (0.20, "#D73027"),
            (0.30, "#F46D43"),
            (0.40, "#FDAE61"),
            (0.50, "#FEE08B"),
            (0.58, "#D9EF8B"),
            (0.66, "#66BD63"),
            (0.74, "#1A9850"),
            (0.82, "#00A6A6"),
            (0.90, "#2C7FB8"),
            (1.00, "#253494"),
        ],
    )

    fig, ax = plt.subplots(
        figsize=(13, 11),
        constrained_layout=True,
    )

    image = ax.imshow(
        matrix,
        cmap=colour_map,
        norm=colour_norm,
        interpolation="nearest",
        aspect="equal",
    )

    ticks = np.arange(len(labels))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.set_xticklabels(
        labels,
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=10,
    )
    ax.set_yticklabels(labels, fontsize=10)

    ax.set_xlabel("Missense variant", labelpad=12)
    ax.set_ylabel("Missense variant", labelpad=12)

    ax.set_title(
        "Pairwise haplotype LD among 11 representative-transcript "
        "missense variants\n"
        "Q9774: 165 cattle / 330 phased chromosome copies; "
        "MAF ≥ 0.01; ARS-UCD1.2",
        fontsize=12,
        pad=18,
    )

    ax.set_xticks(
        np.arange(-0.5, len(labels), 1),
        minor=True,
    )
    ax.set_yticks(
        np.arange(-0.5, len(labels), 1),
        minor=True,
    )
    ax.grid(
        which="minor",
        color="white",
        linestyle="-",
        linewidth=1.2,
    )
    ax.tick_params(which="minor", bottom=False, left=False)

    for row in range(len(labels)):
        for column in range(len(labels)):
            value = matrix[row, column]
            red, green, blue, _ = colour_map(colour_norm(value))
            luminance = (
                0.2126 * red +
                0.7152 * green +
                0.0722 * blue
            )
            text_colour = "white" if luminance < 0.50 else "#222222"

            ax.text(
                column,
                row,
                f"{value:.3f}",
                ha="center",
                va="center",
                fontsize=7.6,
                color=text_colour,
            )

    colour_bar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.046,
        pad=0.04,
        ticks=[
            0.000,
            0.001,
            0.003,
            0.006,
            0.007,
            0.010,
            0.013,
            0.014,
            0.025,
            0.050,
            0.100,
            0.200,
            0.500,
            0.750,
            1.000,
        ],
    )
    colour_bar.set_label(
        "Pairwise LD ($r^2$; piecewise nonlinear colour scale)",
        rotation=270,
        labelpad=28,
    )

    fig.suptitle(
        "LCORL missense linkage disequilibrium",
        fontsize=17,
        fontweight="bold",
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


def main():
    for path in [VCF, EVIDENCE, OLD_MATRIX]:
        if not path.is_file():
            raise FileNotFoundError(f"Required input missing: {path}")

    variants = read_variants()
    haplotypes, sample_count = read_phased_alt_alleles(variants)

    chromosome_copies = 2 * sample_count
    if sample_count != 165 or chromosome_copies != 330:
        raise ValueError(
            f"Expected 165 samples/330 chromosome copies; found "
            f"{sample_count}/{chromosome_copies}"
        )

    matrix = calculate_matrix(variants, haplotypes)

    if not np.allclose(matrix, matrix.T, atol=1e-12):
        raise ValueError("LD matrix is not symmetric")

    if not np.allclose(np.diag(matrix), 1.0, atol=1e-12):
        raise ValueError("LD matrix diagonal is not 1")

    validate_against_old_matrix(variants, matrix)
    write_matrix(variants, matrix)
    plot_heatmap(variants, matrix)

    print("Variants:", len(variants))
    print("Samples:", sample_count)
    print("Phased chromosome copies:", chromosome_copies)
    print("Created:", OUTPUT_MATRIX)
    print("Created:", OUTPUT_PNG)
    print("Created:", OUTPUT_PDF)


if __name__ == "__main__":
    main()
