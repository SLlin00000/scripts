from pathlib import Path
import csv
import math

SUMMARY = Path("summary")
OUT = SUMMARY / "figures"
OUT.mkdir(exist_ok=True)

# ---------- helpers ----------
def svg_start(width, height):
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<style>
text {{ font-family: Arial, sans-serif; fill: black; }}
.title {{ font-size: 22px; font-weight: bold; }}
.label {{ font-size: 15px; }}
.small {{ font-size: 12px; }}
</style>
'''

def svg_end():
    return "</svg>"

# ==========================================================
# Figure 1: Raw vs PASS
# ==========================================================
data = [
    ("SNP Raw", 51784),
    ("SNP PASS", 46081),
    ("Indel Raw", 7130),
    ("Indel PASS", 6818),
]

W, H = 800, 520
left, bottom, top = 100, 430, 80
plot_h = bottom - top
maxv = max(v for _, v in data)

svg = svg_start(W, H)
svg += '<text x="400" y="40" text-anchor="middle" class="title">LCORL variant calling and hard filtering</text>'

bar_w = 110
gap = 55

for i, (name, value) in enumerate(data):
    x = left + i * (bar_w + gap)
    h = value / maxv * plot_h
    y = bottom - h

    svg += f'<rect x="{x}" y="{y}" width="{bar_w}" height="{h}" fill="steelblue"/>'
    svg += f'<text x="{x + bar_w/2}" y="{y - 8}" text-anchor="middle" class="label">{value:,}</text>'
    svg += f'<text x="{x + bar_w/2}" y="{bottom + 25}" text-anchor="middle" class="label">{name}</text>'

svg += f'<line x1="{left-20}" y1="{bottom}" x2="{W-40}" y2="{bottom}" stroke="black"/>'
svg += svg_end()

(OUT / "01_raw_vs_PASS_variants.svg").write_text(svg)


# ==========================================================
# Figure 2: SNP filter categories
# ==========================================================
filter_data = []

with open(SUMMARY / "snp_filter_counts.txt") as f:
    for line in f:
        p = line.strip().split(maxsplit=1)
        if len(p) == 2:
            count = int(p[0])
            name = p[1]
            if name != "PASS":
                filter_data.append((name, count))

filter_data.sort(key=lambda x: x[1], reverse=True)
filter_data = filter_data[:10]

W, H = 900, 600
left = 230
right = 80
top = 90
row_h = 42
maxv = max(v for _, v in filter_data)

svg = svg_start(W, H)
svg += '<text x="450" y="40" text-anchor="middle" class="title">Major SNP hard-filter outcomes</text>'

for i, (name, value) in enumerate(filter_data):
    y = top + i * row_h
    width = value / maxv * (W-left-right)

    svg += f'<text x="{left-10}" y="{y+18}" text-anchor="end" class="label">{name}</text>'
    svg += f'<rect x="{left}" y="{y}" width="{width}" height="25" fill="steelblue"/>'
    svg += f'<text x="{left+width+8}" y="{y+18}" class="label">{value:,}</text>'

svg += svg_end()
(OUT / "02_SNP_filter_categories.svg").write_text(svg)


# ==========================================================
# Figure 3: Allele frequency histogram
# ==========================================================
afs = []

with open(SUMMARY / "snps_PASS_biallelic_AC_AF.tsv") as f:
    for line in f:
        p = line.strip().split()
        if len(p) >= 5:
            try:
                afs.append(float(p[4]))
            except ValueError:
                pass

bins = 40
counts = [0] * bins

for af in afs:
    idx = min(int(af * bins), bins - 1)
    counts[idx] += 1

W, H = 900, 550
left, right, top, bottom = 90, 40, 80, 450
plot_w = W-left-right
plot_h = bottom-top
maxc = max(counts)

svg = svg_start(W, H)
svg += '<text x="450" y="40" text-anchor="middle" class="title">Allele-frequency distribution of PASS biallelic SNPs</text>'

bar_w = plot_w / bins

for i, c in enumerate(counts):
    h = c / maxc * plot_h if maxc else 0
    x = left + i * bar_w
    y = bottom - h
    svg += f'<rect x="{x}" y="{y}" width="{bar_w-1}" height="{h}" fill="steelblue"/>'

svg += f'<line x1="{left}" y1="{bottom}" x2="{W-right}" y2="{bottom}" stroke="black"/>'
svg += f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="black"/>'

for val in [0, 0.25, 0.5, 0.75, 1.0]:
    x = left + val * plot_w
    svg += f'<text x="{x}" y="{bottom+25}" text-anchor="middle" class="small">{val}</text>'

svg += f'<text x="{W/2}" y="{bottom+55}" text-anchor="middle" class="label">Alternate allele frequency</text>'
svg += '<text x="20" y="270" transform="rotate(-90 20 270)" text-anchor="middle" class="label">Number of SNPs</text>'

svg += svg_end()
(OUT / "03_PASS_SNP_allele_frequency.svg").write_text(svg)


# ==========================================================
# Figure 4: SNP density across region
# ==========================================================
positions = []

with open(SUMMARY / "snps_PASS_positions.tsv") as f:
    for line in f:
        p = line.strip().split()
        if len(p) >= 2:
            try:
                positions.append(int(p[1]))
            except ValueError:
                pass

start = 36377806
end = 38557145
bin_size = 100000
nbins = math.ceil((end-start)/bin_size)
density = [0] * nbins

for pos in positions:
    if start <= pos <= end:
        idx = min((pos-start)//bin_size, nbins-1)
        density[idx] += 1

W, H = 950, 550
left, right, top, bottom = 100, 40, 80, 450
plot_w = W-left-right
plot_h = bottom-top
maxd = max(density)

svg = svg_start(W, H)
svg += '<text x="475" y="40" text-anchor="middle" class="title">PASS SNP density across the LCORL region</text>'

points = []

for i, d in enumerate(density):
    x = left + (i/(nbins-1)) * plot_w
    y = bottom - (d/maxd * plot_h if maxd else 0)
    points.append(f"{x},{y}")

svg += f'<polyline points="{" ".join(points)}" fill="none" stroke="steelblue" stroke-width="3"/>'
svg += f'<line x1="{left}" y1="{bottom}" x2="{W-right}" y2="{bottom}" stroke="black"/>'
svg += f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="black"/>'

svg += f'<text x="{left}" y="{bottom+25}" text-anchor="middle" class="small">{start/1e6:.2f} Mb</text>'
svg += f'<text x="{W-right}" y="{bottom+25}" text-anchor="middle" class="small">{end/1e6:.2f} Mb</text>'
svg += f'<text x="{W/2}" y="{bottom+55}" text-anchor="middle" class="label">Chromosome 6 position</text>'
svg += '<text x="20" y="270" transform="rotate(-90 20 270)" text-anchor="middle" class="label">PASS SNPs per 100 kb</text>'

svg += svg_end()
(OUT / "04_PASS_SNP_density_LCORL.svg").write_text(svg)

print("DONE")
print("Figures are in:", OUT)
