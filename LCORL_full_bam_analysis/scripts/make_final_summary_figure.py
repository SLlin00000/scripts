import math
import os

# ============================================================
# INPUTS
# ============================================================

maf_file = "analysis/PASS_biallelic_SNP_AF_MAF.tsv"
density_file = "analysis/PASS_SNP_density_bins.tsv"
gene_file = "analysis/LCORL_gene_span.tsv"
ld_file = "haplotype/LCORL_gene_LD_decay_10kb.tsv"

out_file = "summary/figures/07_LCORL_final_summary.svg"

# ============================================================
# HELPERS
# ============================================================

def read_table(path):
    rows = []
    with open(path) as f:
        header = f.readline().strip().split()
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) != len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return rows

def esc(s):
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )

def text(x, y, s, size=24, anchor="start", weight="normal",
         fill="#111111"):
    return (
        f'<text x="{x}" y="{y}" '
        f'font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" '
        f'text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>'
    )

def line(x1, y1, x2, y2, stroke="#222222", width=2, dash=None):
    extra = ""
    if dash:
        extra = f' stroke-dasharray="{dash}"'
    return (
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
        f'stroke="{stroke}" stroke-width="{width}"{extra}/>'
    )

def rect(x, y, w, h, fill="none", stroke="none", width=1, rx=0):
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{width}" '
        f'rx="{rx}"/>'
    )

def circle(cx, cy, r, fill="#333333"):
    return f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}"/>'

# ============================================================
# READ DATA
# ============================================================

# ---------- MAF ----------
maf_rows = read_table(maf_file)

maf_values = []
for r in maf_rows:
    try:
        v = float(r["MAF"])
        if 0 <= v <= 0.5:
            maf_values.append(v)
    except:
        pass

# 0.05-width bins
maf_edges = [i / 20 for i in range(11)]   # 0, .05 ... .50
maf_counts = [0] * 10

for v in maf_values:
    idx = int(v / 0.05)
    if idx >= 10:
        idx = 9
    maf_counts[idx] += 1

# ---------- density ----------
density_rows = read_table(density_file)

density = []
for r in density_rows:
    try:
        density.append(
            (
                int(r["bin_start"]),
                int(r["bin_end"]),
                int(r["count"])
            )
        )
    except:
        pass

# ---------- gene ----------
gene_rows = read_table(gene_file)

gene_start = int(gene_rows[0]["START"])
gene_end = int(gene_rows[0]["END"])
gene_name = gene_rows[0]["GENE"]

# ---------- LD ----------
ld_rows = read_table(ld_file)

ld = []
for r in ld_rows:
    try:
        b1 = int(r["BIN_START"])
        b2 = int(r["BIN_END"])
        n = int(r["N_PAIRS"])
        r2 = float(r["MEAN_R2"])
        midpoint_kb = ((b1 + b2) / 2) / 1000
        ld.append((midpoint_kb, n, r2))
    except:
        pass

# ============================================================
# CANVAS
# ============================================================

W = 1800
H = 1500

svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" '
    f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
    rect(0, 0, W, H, fill="white")
]

# ============================================================
# MAIN TITLE
# ============================================================

svg.append(
    text(
        W / 2, 60,
        "LCORL variant and haplotype analysis summary",
        size=34,
        anchor="middle",
        weight="bold"
    )
)

svg.append(
    text(
        W / 2, 98,
        "Chromosome 6 | tropical cattle whole-genome sequencing dataset",
        size=19,
        anchor="middle",
        fill="#555555"
    )
)

# ============================================================
# WORKFLOW SUMMARY BOX
# ============================================================

box_x = 100
box_y = 130
box_w = 1600
box_h = 110

svg.append(
    rect(box_x, box_y, box_w, box_h,
         fill="#F7F7F7", stroke="#CCCCCC", width=2, rx=14)
)

workflow = (
    "391 BAMs  →  174 biological samples  →  "
    "46,081 PASS SNPs  →  genotype QC  →  "
    "165 samples / 43,902 SNPs  →  Beagle phasing  →  "
    "2,213 SNPs within LCORL"
)

svg.append(
    text(
        W / 2, 180,
        workflow,
        size=22,
        anchor="middle",
        weight="bold"
    )
)

svg.append(
    text(
        W / 2, 215,
        "QC thresholds: sample missingness ≤ 5%; variant missingness ≤ 10%",
        size=17,
        anchor="middle",
        fill="#555555"
    )
)

# ============================================================
# PANEL CONSTANTS
# ============================================================

panel_x = 100
panel_w = 1600
panel_h = 330

plot_left = panel_x + 115
plot_right = panel_x + panel_w - 65
plot_w = plot_right - plot_left

# ============================================================
# PANEL A — MAF DISTRIBUTION
# ============================================================

py = 285

svg.append(
    rect(panel_x, py, panel_w, panel_h,
         fill="#FFFFFF", stroke="#DDDDDD", width=2, rx=12)
)

svg.append(
    text(panel_x + 25, py + 38, "A", size=28, weight="bold")
)

svg.append(
    text(
        panel_x + 75, py + 38,
        "Minor allele frequency distribution of PASS biallelic SNPs",
        size=24,
        weight="bold"
    )
)

x0 = plot_left
y0 = py + 270
plot_h = 180

svg.append(line(x0, y0, plot_right, y0))
svg.append(line(x0, y0 - plot_h, x0, y0))

max_maf = max(maf_counts) if maf_counts else 1
bar_w = plot_w / len(maf_counts)

for i, c in enumerate(maf_counts):
    bh = (c / max_maf) * plot_h
    bx = x0 + i * bar_w + 8
    bw = bar_w - 16

    svg.append(
        rect(bx, y0 - bh, bw, bh,
             fill="#5B8FF9", stroke="#FFFFFF")
    )

    lab = f"{maf_edges[i]:.2f}–{maf_edges[i+1]:.2f}"
    svg.append(
        text(
            bx + bw / 2, y0 + 28,
            lab,
            size=14,
            anchor="middle"
        )
    )

# y ticks
for k in range(5):
    val = max_maf * k / 4
    yy = y0 - plot_h * k / 4
    svg.append(line(x0 - 6, yy, x0, yy))
    svg.append(
        text(x0 - 12, yy + 5, f"{int(val):,}",
             size=14, anchor="end")
    )

svg.append(
    text(
        (x0 + plot_right) / 2, py + 325,
        "Minor allele frequency (MAF)",
        size=17,
        anchor="middle"
    )
)

svg.append(
    text(
        x0 - 75, py + 175,
        "SNP count",
        size=17,
        anchor="middle"
    )
)

# ============================================================
# PANEL B — SNP DENSITY
# ============================================================

py = 650

svg.append(
    rect(panel_x, py, panel_w, panel_h,
         fill="#FFFFFF", stroke="#DDDDDD", width=2, rx=12)
)

svg.append(
    text(panel_x + 25, py + 38, "B", size=28, weight="bold")
)

svg.append(
    text(
        panel_x + 75, py + 38,
        "PASS SNP density across the LCORL ±1 Mb region",
        size=24,
        weight="bold"
    )
)

x0 = plot_left
y0 = py + 245
plot_h = 145

region_start = density[0][0]
region_end = density[-1][1]

max_density = max(d[2] for d in density)

density_bar_w = plot_w / len(density)

for i, (bs, be, count) in enumerate(density):
    bh = (count / max_density) * plot_h

    svg.append(
        rect(
            x0 + i * density_bar_w + 2,
            y0 - bh,
            density_bar_w - 4,
            bh,
            fill="#61DDAA"
        )
    )

svg.append(line(x0, y0, plot_right, y0))
svg.append(line(x0, y0 - plot_h, x0, y0))

# gene position mapped to x
def density_x(pos):
    return x0 + (
        (pos - region_start) /
        (region_end - region_start)
    ) * plot_w

gx1 = density_x(gene_start)
gx2 = density_x(gene_end)

gene_y = y0 + 45

svg.append(
    rect(
        gx1, gene_y,
        max(gx2 - gx1, 3), 18,
        fill="#E74C3C"
    )
)

svg.append(
    text(
        (gx1 + gx2) / 2,
        gene_y - 8,
        "LCORL",
        size=16,
        anchor="middle",
        weight="bold",
        fill="#B03A2E"
    )
)

# x labels
positions = [
    region_start,
    gene_start,
    gene_end,
    region_end
]

labels = [
    f"{region_start/1e6:.2f} Mb",
    f"{gene_start/1e6:.2f}",
    f"{gene_end/1e6:.2f}",
    f"{region_end/1e6:.2f} Mb"
]

for p, lab in zip(positions, labels):
    xx = density_x(p)
    svg.append(line(xx, y0, xx, y0 + 7))
    svg.append(
        text(xx, y0 + 28, lab,
             size=14, anchor="middle")
    )

svg.append(
    text(
        (x0 + plot_right) / 2,
        py + 315,
        "Genomic position on chromosome 6",
        size=17,
        anchor="middle"
    )
)

svg.append(
    text(
        x0 - 75, py + 170,
        "SNPs / 100 kb",
        size=17,
        anchor="middle"
    )
)

# ============================================================
# PANEL C — LD DECAY
# ============================================================

py = 1015

svg.append(
    rect(panel_x, py, panel_w, panel_h,
         fill="#FFFFFF", stroke="#DDDDDD", width=2, rx=12)
)

svg.append(
    text(panel_x + 25, py + 38, "C", size=28, weight="bold")
)

svg.append(
    text(
        panel_x + 75, py + 38,
        "LD decay within the phased LCORL gene region",
        size=24,
        weight="bold"
    )
)

x0 = plot_left
y0 = py + 250
plot_h = 155

max_x = max(x[0] for x in ld)
max_r2 = max(x[2] for x in ld)

# Give y some breathing room
ymax = max(0.04, max_r2 * 1.15)

def ld_x(kb):
    return x0 + (kb / max_x) * plot_w

def ld_y(r2):
    return y0 - (r2 / ymax) * plot_h

svg.append(line(x0, y0, plot_right, y0))
svg.append(line(x0, y0 - plot_h, x0, y0))

# horizontal guides
for i in range(5):
    val = ymax * i / 4
    yy = ld_y(val)

    svg.append(
        line(
            x0, yy,
            plot_right, yy,
            stroke="#EEEEEE",
            width=1
        )
    )

    svg.append(
        text(
            x0 - 12, yy + 5,
            f"{val:.3f}",
            size=14,
            anchor="end"
        )
    )

# polyline
points = []

for kb, n, r2 in ld:
    xx = ld_x(kb)
    yy = ld_y(r2)
    points.append(f"{xx:.2f},{yy:.2f}")

svg.append(
    f'<polyline points="{" ".join(points)}" '
    f'fill="none" stroke="#5D7092" stroke-width="4"/>'
)

for kb, n, r2 in ld:
    svg.append(
        circle(
            ld_x(kb),
            ld_y(r2),
            5,
            fill="#5D7092"
        )
    )

# x ticks every ~20 kb
tick_step = 20

tick = 0
while tick <= max_x + 1:
    xx = ld_x(tick)
    svg.append(line(xx, y0, xx, y0 + 7))
    svg.append(
        text(
            xx, y0 + 27,
            f"{int(tick)}",
            size=14,
            anchor="middle"
        )
    )
    tick += tick_step

svg.append(
    text(
        (x0 + plot_right) / 2,
        py + 315,
        "Distance between SNP pairs (kb)",
        size=17,
        anchor="middle"
    )
)

svg.append(
    text(
        x0 - 75, py + 175,
        "Mean r²",
        size=17,
        anchor="middle"
    )
)

# small note in unused upper-right area
svg.append(
    text(
        plot_right - 15,
        py + 75,
        f"Valid LD pairs: 2,359,878",
        size=15,
        anchor="end",
        fill="#555555"
    )
)

svg.append(
    text(
        plot_right - 15,
        py + 100,
        f"10 kb distance bins",
        size=15,
        anchor="end",
        fill="#555555"
    )
)

# ============================================================
# FOOTNOTE
# ============================================================

svg.append(
    text(
        W / 2, 1430,
        "LCORL coordinates: chr6:37,377,806–37,556,970 (ARS-UCD1.2)",
        size=16,
        anchor="middle",
        fill="#555555"
    )
)

svg.append(
    text(
        W / 2, 1460,
        "Phasing: Beagle 5.4 | Pairwise haplotype LD: VCFtools 0.1.16",
        size=16,
        anchor="middle",
        fill="#555555"
    )
)

svg.append("</svg>")

# ============================================================
# WRITE
# ============================================================

os.makedirs(os.path.dirname(out_file), exist_ok=True)

with open(out_file, "w") as f:
    f.write("\n".join(svg))

print("DONE")
print("Final summary figure:")
print(out_file)

