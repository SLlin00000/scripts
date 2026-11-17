import math

# ---------- input files ----------
pos_file = "summary/snps_PASS_positions.tsv"
gene_file = "analysis/LCORL_gene_span.tsv"

# ---------- output files ----------
out_svg = "summary/figures/05_PASS_SNP_density_with_LCORL.svg"
out_tsv = "analysis/PASS_SNP_density_bins.tsv"

# ---------- read gene span ----------
with open(gene_file) as f:
    header = next(f)
    row = next(f).strip().split("\t")
    chrom = row[0]
    gene_start = int(row[1])
    gene_end = int(row[2])
    gene_name = row[3]

# define plotting region = gene ± 1 Mb
region_start = gene_start - 1000000
region_end = gene_end + 1000000
region_len = region_end - region_start + 1

# ---------- read SNP positions ----------
positions = []
with open(pos_file) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        c = parts[0]
        p = int(parts[1])
        if c == chrom and region_start <= p <= region_end:
            positions.append(p)

# ---------- binning ----------
bin_size = 100000   # 100 kb
n_bins = math.ceil(region_len / bin_size)
counts = [0] * n_bins

for p in positions:
    idx = (p - region_start) // bin_size
    if 0 <= idx < n_bins:
        counts[idx] += 1

# ---------- save bin table ----------
with open(out_tsv, "w") as out:
    out.write("bin_start\tbin_end\tcount\n")
    for i, c in enumerate(counts):
        start = region_start + i * bin_size
        end = min(start + bin_size - 1, region_end)
        out.write(f"{start}\t{end}\t{c}\n")

# ---------- draw simple SVG ----------
W = 1200
H = 700
left = 90
right = 40
top = 70
bottom = 120

plot_w = W - left - right
plot_h = H - top - bottom

max_count = max(counts) if max(counts) > 0 else 1
bar_w = plot_w / n_bins

def x_scale(pos):
    return left + (pos - region_start) / (region_end - region_start) * plot_w

def y_scale(count):
    return top + plot_h - (count / max_count) * plot_h

svg = []
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}">')
svg.append(f'<rect width="100%" height="100%" fill="white"/>')

# title
svg.append(f'<text x="{W/2}" y="35" text-anchor="middle" font-size="22" font-family="Arial">PASS SNP density across LCORL ±1 Mb region</text>')
svg.append(f'<text x="{W/2}" y="58" text-anchor="middle" font-size="14" font-family="Arial">Chromosome {chrom} | total PASS SNPs in region = {len(positions)}</text>')

# axes
svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="black" stroke-width="2"/>')
svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="black" stroke-width="2"/>')

# y ticks
for i in range(6):
    val = round(max_count * i / 5)
    y = y_scale(val)
    svg.append(f'<line x1="{left-5}" y1="{y}" x2="{left}" y2="{y}" stroke="black"/>')
    svg.append(f'<text x="{left-10}" y="{y+5}" text-anchor="end" font-size="12" font-family="Arial">{val}</text>')

# bars
for i, c in enumerate(counts):
    x = left + i * bar_w
    y = y_scale(c)
    h = top + plot_h - y
    svg.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w-2:.2f}" height="{h:.2f}" fill="#5B8FF9" stroke="none"/>')

# gene span
gx1 = x_scale(gene_start)
gx2 = x_scale(gene_end)
gene_y = top + plot_h + 25
svg.append(f'<rect x="{gx1:.2f}" y="{gene_y}" width="{gx2-gx1:.2f}" height="18" fill="#E74C3C"/>')
svg.append(f'<text x="{(gx1+gx2)/2:.2f}" y="{gene_y-5}" text-anchor="middle" font-size="13" font-family="Arial">{gene_name}</text>')

# x ticks
tick_vals = [region_start, gene_start, gene_end, region_end]
tick_labels = [str(region_start), f"{gene_start}\n(gene start)", f"{gene_end}\n(gene end)", str(region_end)]

for xv, lab in zip(tick_vals, tick_labels):
    x = x_scale(xv)
    svg.append(f'<line x1="{x:.2f}" y1="{top+plot_h}" x2="{x:.2f}" y2="{top+plot_h+6}" stroke="black"/>')
    if "\n" in lab:
        a, b = lab.split("\n")
        svg.append(f'<text x="{x:.2f}" y="{top+plot_h+22}" text-anchor="middle" font-size="12" font-family="Arial">{a}</text>')
        svg.append(f'<text x="{x:.2f}" y="{top+plot_h+38}" text-anchor="middle" font-size="11" font-family="Arial">{b}</text>')
    else:
        svg.append(f'<text x="{x:.2f}" y="{top+plot_h+22}" text-anchor="middle" font-size="12" font-family="Arial">{lab}</text>')

# axis labels
svg.append(f'<text x="{W/2}" y="{H-25}" text-anchor="middle" font-size="16" font-family="Arial">Genomic position (bp)</text>')
svg.append(f'<text x="25" y="{H/2}" transform="rotate(-90 25,{H/2})" text-anchor="middle" font-size="16" font-family="Arial">PASS SNP count per 100 kb bin</text>')

svg.append('</svg>')

with open(out_svg, "w") as f:
    f.write("\n".join(svg))

print("DONE")
print("SVG:", out_svg)
print("TSV:", out_tsv)
print("PASS SNPs used:", len(positions))
print("Bins:", n_bins)

