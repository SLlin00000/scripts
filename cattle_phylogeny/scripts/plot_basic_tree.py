from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tree_file = "results/cattle1217_LCORL_NJ_tree.newick"

tree = Phylo.read(tree_file, "newick")

print("Tips:", len(tree.get_terminals()))
print("Internal nodes:", len(tree.get_nonterminals()))

# 为了显示更整齐一点
tree.ladderize()

# 1217 tips 很多，所以先不显示标签
fig = plt.figure(figsize=(14, 40), dpi=300)
ax = fig.add_subplot(1, 1, 1)

Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    label_func=lambda x: None
)

ax.set_title("LCORL NJ tree (1217 cattle, unlabeled overview)", fontsize=14)

plt.tight_layout()
plt.savefig("results/cattle1217_LCORL_NJ_tree_basic_unlabeled.png", dpi=300, bbox_inches="tight")
plt.savefig("results/cattle1217_LCORL_NJ_tree_basic_unlabeled.pdf", bbox_inches="tight")

print("Written:")
print("  results/cattle1217_LCORL_NJ_tree_basic_unlabeled.png")
print("  results/cattle1217_LCORL_NJ_tree_basic_unlabeled.pdf")
