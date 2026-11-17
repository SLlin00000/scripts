from Bio import Phylo
import csv
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

TREE = "results/cattle1217_LCORL_NJ_tree.newick"
META = "metadata/tree1217_annotation_clean.tsv"

OUT_PNG = "results/cattle1217_LCORL_NJ_tree_integrated_v4_rightcolumn.png"
OUT_PDF = "results/cattle1217_LCORL_NJ_tree_integrated_v4_rightcolumn.pdf"

# -----------------------------
# Read metadata
# -----------------------------
meta = {}

with open(META) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        meta[row["sample_id"]] = {
            "source": row["source"],
            "genetic_ancestry": row["genetic_ancestry"],
            "production_purpose": row["production_purpose"],
            "breed_composition": row["breed_composition"]
        }

# -----------------------------
# Colour palettes
# -----------------------------
SOURCE_COLORS = {
    "LOCAL174": "#1B7F3B",
    "RUN8": "#BFC7CE"
}

ANCESTRY_COLORS = {
    "Taurine": "#3B82C4",
    "Indicine": "#E68613",
    "Admixed": "#8E63CE",
    "Needs manual review": "#C99BD8",
    "Metadata not accessible": "#D9D9D9"
}

PURPOSE_COLORS = {
    "Beef": "#C44E52",
    "Dairy": "#4DB6AC",
    "Beef + dairy": "#DDAA33",
    "Purpose unsure": "#8C8C8C",
    "Metadata not accessible": "#D9D9D9"
}

COMPOSITION_COLORS = {
    "Not flagged as cross": "#7EA6E0",
    "Crossbred": "#7F7F7F",
    "Dairy cross": "#5AAE61",
    "Local breed": "#CC79A7",
    "Metadata not accessible": "#D9D9D9"
}

# -----------------------------
# Read tree
# -----------------------------
tree = Phylo.read(TREE, "newick")

try:
    tree.root_at_midpoint()
except Exception:
    pass

tree.ladderize()

terminals = tree.get_terminals()
n = len(terminals)

print("Tree tips:", n)

# -----------------------------
# Tip angles
# -----------------------------
tip_theta = {}

for i, tip in enumerate(terminals):
    tip_theta[tip] = 2.0 * math.pi * i / n

def circ_mean(angles):
    if not angles:
        return 0.0
    x = sum(math.cos(a) for a in angles)
    y = sum(math.sin(a) for a in angles)
    return math.atan2(y, x) % (2.0 * math.pi)

node_theta = {}

for clade in tree.find_clades(order="postorder"):
    if clade.is_terminal():
        node_theta[clade] = tip_theta[clade]
    else:
        node_theta[clade] = circ_mean(
            [node_theta[c] for c in clade.clades]
        )

# -----------------------------
# Radial positions
# -----------------------------
depths = tree.depths()

if not depths or max(depths.values()) == 0:
    depths = tree.depths(unit_branch_lengths=True)

max_depth = max(depths.values())

node_r = {
    clade: depths.get(clade, 0.0) / max_depth
    for clade in tree.find_clades()
}

# -----------------------------
# Figure layout
# -----------------------------
fig = plt.figure(figsize=(22, 16), dpi=300)

# Tree occupies left ~67%; right column is reserved entirely for guide + legends
ax = fig.add_axes([0.03, 0.08, 0.64, 0.84], projection="polar")

ax.set_theta_direction(-1)
ax.set_theta_offset(math.pi / 2.0)
ax.set_axis_off()

# -----------------------------
# Draw tree branches
# -----------------------------
for parent in tree.find_clades(order="level"):
    for child in parent.clades:
        rp = node_r[parent]
        rc = node_r[child]
        tp = node_theta[parent]
        tc = node_theta[child]

        ax.plot(
            [tc, tc],
            [rp, rc],
            linewidth=0.35,
            color="#5A5A5A",
            alpha=0.85,
            zorder=1
        )

        dtheta = (tc - tp + math.pi) % (2.0 * math.pi) - math.pi
        steps = 25

        theta_arc = [
            tp + dtheta * j / steps
            for j in range(steps + 1)
        ]

        r_arc = [rp] * (steps + 1)

        ax.plot(
            theta_arc,
            r_arc,
            linewidth=0.35,
            color="#5A5A5A",
            alpha=0.85,
            zorder=1
        )

# -----------------------------
# Annotation rings
# -----------------------------
width = 2.0 * math.pi / n

ring1_bottom = 1.035
ring2_bottom = 1.085
ring3_bottom = 1.135
ring4_bottom = 1.185
ring_height = 0.036

for tip in terminals:
    sid = tip.name
    theta = tip_theta[tip]

    m = meta.get(
        sid,
        {
            "source": "RUN8",
            "genetic_ancestry": "Metadata not accessible",
            "production_purpose": "Metadata not accessible",
            "breed_composition": "Metadata not accessible"
        }
    )

    ax.bar(
        theta, ring_height,
        width=width * 1.02,
        bottom=ring1_bottom,
        color=SOURCE_COLORS.get(m["source"], "#E5E5E5"),
        edgecolor="none",
        align="edge",
        zorder=2
    )

    ax.bar(
        theta, ring_height,
        width=width * 1.02,
        bottom=ring2_bottom,
        color=ANCESTRY_COLORS.get(
            m["genetic_ancestry"], "#D9D9D9"
        ),
        edgecolor="none",
        align="edge",
        zorder=2
    )

    ax.bar(
        theta, ring_height,
        width=width * 1.02,
        bottom=ring3_bottom,
        color=PURPOSE_COLORS.get(
            m["production_purpose"], "#D9D9D9"
        ),
        edgecolor="none",
        align="edge",
        zorder=2
    )

    ax.bar(
        theta, ring_height,
        width=width * 1.02,
        bottom=ring4_bottom,
        color=COMPOSITION_COLORS.get(
            m["breed_composition"], "#D9D9D9"
        ),
        edgecolor="none",
        align="edge",
        zorder=2
    )

ax.set_ylim(0, 1.25)

# -----------------------------
# Right-side guide + legends
# -----------------------------
# Entire right column reserved for text/legend only

# Ring order guide at very top
fig.text(
    0.72,
    0.915,
    "Ring order",
    ha="left",
    va="top",
    fontsize=13,
    fontweight="bold"
)

fig.text(
    0.72,
    0.875,
    "Outer ring\n"
    "↓\n"
    "Inner ring",
    ha="left",
    va="top",
    fontsize=11,
    linespacing=1.5
)

fig.text(
    0.79,
    0.875,
    "1. Breed composition\n"
    "2. Production purpose\n"
    "3. Genetic ancestry\n"
    "4. Source",
    ha="left",
    va="top",
    fontsize=11,
    linespacing=1.5
)

# Legend handles
source_handles = [
    Patch(facecolor=SOURCE_COLORS[k], label=k)
    for k in ["LOCAL174", "RUN8"]
]

ancestry_handles = [
    Patch(facecolor=ANCESTRY_COLORS[k], label=k)
    for k in [
        "Taurine",
        "Indicine",
        "Admixed",
        "Needs manual review",
        "Metadata not accessible"
    ]
]

purpose_handles = [
    Patch(facecolor=PURPOSE_COLORS[k], label=k)
    for k in [
        "Beef",
        "Dairy",
        "Beef + dairy",
        "Purpose unsure",
        "Metadata not accessible"
    ]
]

composition_handles = [
    Patch(facecolor=COMPOSITION_COLORS[k], label=k)
    for k in [
        "Not flagged as cross",
        "Crossbred",
        "Dairy cross",
        "Local breed",
        "Metadata not accessible"
    ]
]

# Sample source
leg1 = fig.legend(
    handles=source_handles,
    title="Sample source",
    loc="upper left",
    bbox_to_anchor=(0.72, 0.69),
    frameon=False,
    fontsize=11,
    title_fontsize=12,
    labelspacing=0.8,
    handlelength=1.5
)
fig.add_artist(leg1)

# Genetic ancestry
leg2 = fig.legend(
    handles=ancestry_handles,
    title="Genetic ancestry",
    loc="upper left",
    bbox_to_anchor=(0.72, 0.52),
    frameon=False,
    fontsize=11,
    title_fontsize=12,
    labelspacing=0.8,
    handlelength=1.5
)
fig.add_artist(leg2)

# Production purpose
leg3 = fig.legend(
    handles=purpose_handles,
    title="Production purpose",
    loc="upper left",
    bbox_to_anchor=(0.72, 0.31),
    frameon=False,
    fontsize=11,
    title_fontsize=12,
    labelspacing=0.8,
    handlelength=1.5
)
fig.add_artist(leg3)

# Breed composition
leg4 = fig.legend(
    handles=composition_handles,
    title="Breed composition",
    loc="upper left",
    bbox_to_anchor=(0.72, 0.08),
    frameon=False,
    fontsize=11,
    title_fontsize=12,
    labelspacing=0.8,
    handlelength=1.5
)
fig.add_artist(leg4)

# -----------------------------
# Title / footnotes
# -----------------------------
fig.suptitle(
    "LCORL-region genetic structure in 1,217 cattle",
    fontsize=18,
    y=0.975
)

fig.text(
    0.34,
    0.035,
    "Neighbour-joining tree based on 4,056 LD-pruned SNPs "
    "(MAF ≥ 0.01; genotype missingness < 5%)",
    ha="center",
    fontsize=11
)

fig.text(
    0.34,
    0.017,
    "LOCAL174 annotations for ancestry, production purpose, and breed composition "
    "were not accessible at the time of analysis.",
    ha="center",
    fontsize=10
)

plt.savefig(
    OUT_PNG,
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    OUT_PDF,
    bbox_inches="tight"
)

print("Written:")
print(OUT_PNG)
print(OUT_PDF)
