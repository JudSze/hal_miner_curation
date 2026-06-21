import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("model_evaluation/hit_lineages/dimetal_lineage_with_counts.tsv", sep="\t")

phylum_counts = df.groupby("phylum")["n_sequences"].sum().sort_values(ascending=False)

fig, ax = plt.subplots(figsize=(8, 24))
bars = ax.barh(phylum_counts.index[::-1], phylum_counts.values[::-1],
               color="#4878CF", edgecolor="white")

for bar, val in zip(bars, phylum_counts.values[::-1]):
    ax.text(bar.get_width() + 2, bar.get_y() + bar.get_height() / 2,
            str(int(val)), va="center", fontsize=8)

ax.set_xlabel("Number of sequences")
ax.set_title("Conventional FDH hit sequences — distribution by phylum\n"
             f"Total: {phylum_counts.sum()}", fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(axis='y', labelsize=8)
plt.tight_layout()
plt.savefig("fdh_conventional.png", bbox_inches="tight", dpi=150)