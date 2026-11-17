import os

INDIR = "haplotype/segment_haplotypes"
OUT = "summary/figures/08_LCORL_candidate_haplotype_frequencies.svg"

segments = ["S1", "S6", "S3", "S4", "S5"]

data = {}

for seg in segments:
    path = f"{INDIR}/{seg}_haplotype_frequencies.tsv"

    rows = []

    with open(path) as f:
        header = f.readline()

        for line in f:
            if not line.strip():
                continue

            p = line.rstrip().split("\t")

            hap = p[2]
            count = int(p[3])
            freq = float(p[4])

            rows.append((hap, count, freq))

    data[seg] = rows


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


W = 1600
H = 1200

svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" '
    f'width="{W}" height="{H}" '
    f'viewBox="0 0 {W} {H}">',
    rect(0, 0, W, H, fill="white")
]

svg.append(
    text(
        W/2, 55,
        "Candidate LCORL haplotype frequencies",
        size=32,
        anchor="middle",
        weight="bold"
    )
)

svg.append(
    text(
        W/2, 90,
        "Phased common SNPs (MAF ≥ 0.05); 330 chromosomes",
        size=18,
        anchor="middle",
        fill="#555555"
    )
)

panel_top = 135
panel_h = 185
left = 120
bar_start = 420
bar_max_w = 900

for idx, seg in enumerate(segments):

    y = panel_top + idx * 205

    svg.append(
        rect(
            80, y-30, 1440, panel_h,
            fill="#FAFAFA",
            stroke="#DDDDDD",
            sw=1.5,
            rx=10
        )
    )

    priority = (
        "Primary candidate"
        if seg in ["S1", "S6"]
        else "Secondary candidate"
    )

    svg.append(
        text(
            left, y,
            f"{seg} — {priority}",
            size=22,
            weight="bold"
        )
    )

    rows = data[seg]

    # show common haplotypes individually;
    # pool rare haplotypes into Other
    common = []
    other = 0.0

    for hap, count, freq in rows:
        if freq >= 0.05:
            common.append((hap, freq))
        else:
            other += freq

    if other > 0:
        common.append(("Other", other))

    row_y = y + 38

    for hap, freq in common:

        svg.append(
            text(
                left,
                row_y + 18,
                hap,
                size=17
            )
        )

        bw = freq * bar_max_w

        svg.append(
            rect(
                bar_start,
                row_y,
                bw,
                26,
                fill="#5B8FF9"
            )
        )

        svg.append(
            text(
                bar_start + bw + 12,
                row_y + 19,
                f"{freq*100:.1f}%",
                size=16
            )
        )

        row_y += 38


svg.append(
    text(
        W/2,
        H-35,
        "Primary candidates were prioritised by haplotype diversity and frequency balance.",
        size=16,
        anchor="middle",
        fill="#555555"
    )
)

svg.append("</svg>")

os.makedirs(os.path.dirname(OUT), exist_ok=True)

with open(OUT, "w") as f:
    f.write("\n".join(svg))

print("DONE")
print(OUT)

