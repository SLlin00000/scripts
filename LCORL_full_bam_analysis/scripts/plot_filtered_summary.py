import os
import math

OUT = "summary/figures/09_LCORL_filtered_summary.svg"

# --------------------------------------------------
# Known results
# --------------------------------------------------

raw_snps = 51784
pass_snps = 46081
filtered_total = 5703

filter_counts = [
    ("MQ < 40", 2661),
    ("SOR > 3", 2277),
    ("QD < 2", 1776),
    ("FS > 60", 239),
]

pass_pct = pass_snps / raw_snps * 100
filtered_pct = filtered_total / raw_snps * 100

# Other confirmed project results
qc_snps = 43902
common_gene_snps = 720

# --------------------------------------------------
# SVG helpers
# --------------------------------------------------

def esc(x):
    return (
        str(x)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )

def text(x, y, s, size=18, anchor="start",
         weight="normal", fill="#111111"):
    return (
        f'<text x="{x}" y="{y}" '
        f'font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" '
        f'font-weight="{weight}" '
        f'text-anchor="{anchor}" '
        f'fill="{fill}">{esc(s)}</text>'
    )

def rect(x, y, w, h, fill="none",
         stroke="none", sw=1, rx=0):
    return (
        f'<rect x="{x}" y="{y}" '
        f'width="{w}" height="{h}" '
        f'fill="{fill}" stroke="{stroke}" '
        f'stroke-width="{sw}" rx="{rx}"/>'
    )

def line(x1, y1, x2, y2, stroke="#333333", sw=2):
    return (
        f'<line x1="{x1}" y1="{y1}" '
        f'x2="{x2}" y2="{y2}" '
        f'stroke="{stroke}" stroke-width="{sw}"/>'
    )

W = 1700
H = 1100

svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" '
    f'width="{W}" height="{H}" '
    f'viewBox="0 0 {W} {H}">',
    rect(0, 0, W, H, fill="white")
]

# --------------------------------------------------
# Title
# --------------------------------------------------

svg.append(
    text(
        W/2, 55,
        "LCORL variant filtering and downstream SNP selection",
        size=32,
        anchor="middle",
        weight="bold"
    )
)

svg.append(
    text(
        W/2, 90,
        "Full BAM cohort; Chr6 LCORL region",
        size=18,
        anchor="middle",
        fill="#555555"
    )
)

# --------------------------------------------------
# LEFT PANEL: filtering workflow
# --------------------------------------------------

panel_x = 70
panel_y = 135
panel_w = 690
panel_h = 880

svg.append(
    rect(
        panel_x, panel_y, panel_w, panel_h,
        fill="#FAFAFA",
        stroke="#DDDDDD",
        sw=1.5,
        rx=12
    )
)

svg.append(
    text(
        panel_x + 30,
        panel_y + 45,
        "A. Variant filtering workflow",
        size=24,
        weight="bold"
    )
)

box_x = panel_x + 145
box_w = 400
box_h = 95

# Raw SNPs
y1 = panel_y + 110
svg.append(rect(box_x, y1, box_w, box_h,
                fill="#F3F6FA", stroke="#B8C2CC", sw=1.5, rx=10))
svg.append(text(box_x + box_w/2, y1 + 38,
                f"{raw_snps:,} raw SNPs",
                size=26, anchor="middle", weight="bold"))
svg.append(text(box_x + box_w/2, y1 + 68,
                "Before SNP hard filtering",
                size=16, anchor="middle", fill="#555555"))

# Arrow
svg.append(line(box_x + box_w/2, y1 + box_h,
                box_x + box_w/2, y1 + box_h + 55))
svg.append(text(box_x + box_w/2 + 20, y1 + box_h + 35,
                "GATK hard filtering",
                size=16, fill="#555555"))

# PASS
y2 = y1 + 170
svg.append(rect(box_x, y2, box_w, box_h,
                fill="#F3F6FA", stroke="#B8C2CC", sw=1.5, rx=10))
svg.append(text(box_x + box_w/2, y2 + 38,
                f"{pass_snps:,} PASS SNPs",
                size=26, anchor="middle", weight="bold"))
svg.append(text(box_x + box_w/2, y2 + 68,
                f"{pass_pct:.1f}% of raw SNPs",
                size=16, anchor="middle", fill="#555555"))

# Filtered offshoot
svg.append(line(box_x + box_w, y2 + box_h/2,
                box_x + box_w + 80, y2 + box_h/2))
svg.append(text(box_x + box_w + 95, y2 + 35,
                f"{filtered_total:,} filtered",
                size=18, weight="bold"))
svg.append(text(box_x + box_w + 95, y2 + 63,
                f"({filtered_pct:.1f}%)",
                size=16, fill="#555555"))

# Arrow
svg.append(line(box_x + box_w/2, y2 + box_h,
                box_x + box_w/2, y2 + box_h + 55))
svg.append(text(box_x + box_w/2 + 20, y2 + box_h + 35,
                "biallelic + missingness QC",
                size=16, fill="#555555"))

# QC SNPs
y3 = y2 + 170
svg.append(rect(box_x, y3, box_w, box_h,
                fill="#F3F6FA", stroke="#B8C2CC", sw=1.5, rx=10))
svg.append(text(box_x + box_w/2, y3 + 38,
                f"{qc_snps:,} QC SNPs",
                size=26, anchor="middle", weight="bold"))
svg.append(text(box_x + box_w/2, y3 + 68,
                "165 samples retained",
                size=16, anchor="middle", fill="#555555"))

# Arrow
svg.append(line(box_x + box_w/2, y3 + box_h,
                box_x + box_w/2, y3 + box_h + 55))
svg.append(text(box_x + box_w/2 + 20, y3 + box_h + 35,
                "LCORL gene body + MAF ≥ 0.05",
                size=16, fill="#555555"))

# Common SNPs
y4 = y3 + 170
svg.append(rect(box_x, y4, box_w, box_h,
                fill="#F3F6FA", stroke="#B8C2CC", sw=1.5, rx=10))
svg.append(text(box_x + box_w/2, y4 + 38,
                f"{common_gene_snps:,} common phased SNPs",
                size=24, anchor="middle", weight="bold"))
svg.append(text(box_x + box_w/2, y4 + 68,
                "Used for local LD / haplotype analysis",
                size=16, anchor="middle", fill="#555555"))

# --------------------------------------------------
# RIGHT PANEL: filter criteria counts
# --------------------------------------------------

rx = 815
ry = 135
rw = 815
rh = 880

svg.append(
    rect(
        rx, ry, rw, rh,
        fill="#FAFAFA",
        stroke="#DDDDDD",
        sw=1.5,
        rx=12
    )
)

svg.append(
    text(
        rx + 30,
        ry + 45,
        "B. Main hard-filter failure criteria",
        size=24,
        weight="bold"
    )
)

svg.append(
    text(
        rx + 30,
        ry + 78,
        "Counts are non-exclusive; one SNP may fail multiple criteria.",
        size=16,
        fill="#555555"
    )
)

max_count = max(v for _, v in filter_counts)

bar_x = rx + 210
bar_max_w = 500
start_y = ry + 145
row_gap = 115

for i, (label, value) in enumerate(filter_counts):
    yy = start_y + i * row_gap

    svg.append(
        text(
            rx + 35,
            yy + 31,
            label,
            size=20,
            weight="bold"
        )
    )

    bw = value / max_count * bar_max_w

    svg.append(
        rect(
            bar_x,
            yy,
            bw,
            42,
            fill="#5B8FF9"
        )
    )

    svg.append(
        text(
            bar_x + bw + 15,
            yy + 29,
            f"{value:,}",
            size=18
        )
    )

# Zero-count criteria note
note_y = start_y + 4 * row_gap + 25

svg.append(
    text(
        rx + 35,
        note_y,
        "Additional configured filters with zero failures:",
        size=18,
        weight="bold"
    )
)

svg.append(
    text(
        rx + 35,
        note_y + 34,
        "MQRankSum < -12.5; ReadPosRankSum < -8; QUAL < 30",
        size=17,
        fill="#555555"
    )
)

# Interpretive note
svg.append(
    rect(
        rx + 35,
        note_y + 85,
        rw - 70,
        145,
        fill="#FFFFFF",
        stroke="#D8D8D8",
        sw=1.2,
        rx=8
    )
)

svg.append(
    text(
        rx + 60,
        note_y + 120,
        "Interpretation",
        size=20,
        weight="bold"
    )
)

svg.append(
    text(
        rx + 60,
        note_y + 153,
        "Most excluded SNPs were affected by mapping quality,",
        size=17
    )
)

svg.append(
    text(
        rx + 60,
        note_y + 181,
        "strand-bias, or quality-by-depth criteria.",
        size=17
    )
)

# footer
svg.append(
    text(
        W/2,
        H - 28,
        "Hard-filter counts are descriptive and reflect the filters already applied to the SNP VCF.",
        size=15,
        anchor="middle",
        fill="#666666"
    )
)

svg.append("</svg>")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

with open(OUT, "w") as f:
    f.write("\n".join(svg))

print("DONE")
print(OUT)

