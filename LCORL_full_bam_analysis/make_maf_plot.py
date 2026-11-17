from pathlib import Path

infile = Path("analysis/PASS_biallelic_SNP_AF_MAF.tsv")
outfile = Path("analysis/PASS_biallelic_SNP_MAF.svg")

mafs = []

with open(infile) as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:
            mafs.append(float(parts[5]))

# 10 bins: 0–0.05, 0.05–0.10 ... 0.45–0.50
bins = 10
counts = [0] * bins

for maf in mafs:
    idx = min(int(maf / 0.05), 9)
    counts[idx] += 1

W, H = 950, 560
left, right, top, bottom = 100, 40, 80, 450
plot_w = W - left - right
plot_h = bottom - top
maxc = max(counts)

svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}">
<rect width="100%" height="100%" fill="white"/>
<style>
text {{ font-family: Arial, sans-serif; fill: black; }}
.title {{ font-size: 22px; font-weight: bold; }}
.label {{ font-size: 15px; }}
.small {{ font-size: 12px; }}
</style>
<text x="475" y="40" text-anchor="middle" class="title">
Minor allele frequency distribution of PASS biallelic SNPs
</text>
'''

bar_w = plot_w / bins

for i, c in enumerate(counts):
    h = c / maxc * plot_h
    x = left + i * bar_w
    y = bottom - h

    svg += f'<rect x="{x}" y="{y}" width="{bar_w-3}" height="{h}" fill="steelblue"/>'
    svg += f'<text x="{x + bar_w/2}" y="{y-6}" text-anchor="middle" class="small">{c}</text>'

    lo = i * 0.05
    hi = (i + 1) * 0.05
    svg += f'<text x="{x + bar_w/2}" y="{bottom+25}" text-anchor="middle" class="small">{lo:.2f}-{hi:.2f}</text>'

svg += f'<line x1="{left}" y1="{bottom}" x2="{W-right}" y2="{bottom}" stroke="black"/>'
svg += f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="black"/>'

svg += f'<text x="{W/2}" y="{bottom+60}" text-anchor="middle" class="label">Minor allele frequency</text>'
svg += '<text x="25" y="270" transform="rotate(-90 25 270)" text-anchor="middle" class="label">Number of SNPs</text>'

svg += '</svg>'

outfile.write_text(svg)

print("DONE")
print("Bin counts:", counts)
print(outfile)

