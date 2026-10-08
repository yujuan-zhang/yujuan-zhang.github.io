"""Step 7: figures for "do the 7 barcodes change the conclusions?" (reads step 4 and step 6 tables).

figures/07_proportion_impact.png  per cell type: largest change in % from (a) the cell set,
                                  (b) the clustering seed, (c) the count matrix
figures/07_cdc_markers.png        cDC markers: % of cDC expressing, Cell Ranger vs Simpleaf+QCatch cell set
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES, TABLES

INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


# --- 1. three sources of change per cell type ---------------------------------------------
st = pd.read_csv(TABLES / "annotation_stability.tsv", sep="\t")
w = st.pivot_table(index=["counts", "seed", "cell_set"], columns="cell_type", values="pct", fill_value=0)
sets = w.xs("Simpleaf+QCatch cells", level="cell_set") - w.xs("Cell Ranger cells", level="cell_set")
cr_set = w.xs("Cell Ranger cells", level="cell_set")
seed = cr_set.groupby(level="counts").agg(lambda s: s.max() - s.min()).max()           # max over the two count matrices
counts = (cr_set.xs("simpleaf_counts", level="counts") - cr_set.xs("cellranger_counts", level="counts")).abs().max()
cell_set = sets.abs().max()
order = cr_set.xs("cellranger_counts", level="counts").loc[0].sort_values().index      # by abundance
d = pd.DataFrame({"cell set": cell_set, "clustering seed": seed, "count matrix": counts}).loc[order]
d.round(3).to_csv(TABLES / "proportion_impact_by_type.tsv", sep="\t")
print(d.round(2))

labels = {"Low quality (no markers)": "Low quality"}
fig, ax = plt.subplots(figsize=(7.2, 4.6))
y = np.arange(len(d))
for k, (col, color) in enumerate([("cell set", BLUE), ("clustering seed", ORANGE), ("count matrix", AQUA)]):
    ax.barh(y + (1 - k) * 0.26, d[col], height=0.22, color=color, label=col, zorder=2)
for yi, v in zip(y, d["cell set"]):
    ax.text(v + 0.12, yi + 0.26, f"{v:.2f}", va="center", fontsize=7.5, color=INK)
ax.set_yticks(y, [labels.get(t, t) for t in d.index], fontsize=9, color=INK)
ax.set_xlabel("Largest change in proportion (percentage points)", color=INK, fontsize=9)
ax.grid(True, axis="x", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=9, loc="lower right", title="Only this changes", title_fontsize=9)
style(ax)
fig.tight_layout(); fig.savefig(FIGURES / "07_proportion_impact.png", dpi=200, facecolor=SURFACE)

# --- 2. cDC markers under the two cell sets -----------------------------------------------------
g = pd.read_csv(TABLES / "marker_conclusion_genes.tsv", sep="\t")
g = g[(g.cell_type == "cDC") & (~g.after_qc)]
genes = ["CST3", "FCER1A", "CLEC10A", "CD1C"]
fig, ax = plt.subplots(figsize=(6.4, 2.6))
for i, gene in enumerate(genes[::-1]):
    a = g[(g.gene == gene) & (g.cell_set == "Cell Ranger")].iloc[0]
    b = g[(g.gene == gene) & (g.cell_set == "Simpleaf+QCatch")].iloc[0]
    ax.plot([a.pct_in_type, b.pct_in_type], [i, i], color=MUTED, linewidth=1.5, zorder=1)
    ax.scatter(a.pct_in_type, i, s=60, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3,
               label="Cell Ranger (16 cDC)" if i == 0 else None)
    ax.scatter(b.pct_in_type, i, s=60, marker="D", color=AQUA, edgecolor=SURFACE, linewidth=1.5, zorder=3,
               label="Simpleaf + QCatch (18 cDC)" if i == 0 else None)
    ax.text(103, i, f"padj {a.padj:.2g} → {b.padj:.2g}", va="center", fontsize=8, color=INK)
ax.set_yticks(range(len(genes)), genes[::-1], fontsize=9, color=INK)
ax.set_xlim(45, 128); ax.set_xticks([50, 60, 70, 80, 90, 100])
ax.set_xlabel("% of cDC expressing the gene", color=INK, fontsize=9)
ax.grid(True, axis="x", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=8, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
style(ax)
fig.tight_layout(); fig.savefig(FIGURES / "07_cdc_markers.png", dpi=200, facecolor=SURFACE)
print("figures written")
