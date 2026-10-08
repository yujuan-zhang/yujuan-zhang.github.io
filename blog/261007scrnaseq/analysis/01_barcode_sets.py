"""Step 1: the three barcode groups and per-barcode metrics from BOTH quantifications.

Outputs
  tables/barcode_sets_summary.tsv   counts of shared / CR_only / SA_only
  tables/union_cells_metrics.tsv    one row per called barcode (either tool), metrics from both raw matrices
  tables/raw_universe_summary.tsv   raw-matrix barcode universes and their UMI distribution
  figures/01_knee_plots.png         barcode-rank plots, disputed barcodes marked
  figures/01_umi_concordance.png    per-barcode UMI, Cell Ranger vs Simpleaf
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES, GROUP_ORDER, TABLES, load_all, load_sa_qcatch_obs, mito_mask, umi_per_bc

COLORS = {"shared": "#2a78d6", "CR_only": "#eb6834", "SA_only": "#1baf7a"}
MARKERS = {"shared": "o", "CR_only": "s", "SA_only": "D"}
LABELS = {"shared": "Both tools", "CR_only": "Cell Ranger only", "SA_only": "Simpleaf + QCatch only"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

cr, sa, cr_cells, sa_cells = load_all()

# --- 1. groups ---------------------------------------------------------------
union = sorted(cr_cells | sa_cells)
group = pd.Series(
    ["shared" if b in cr_cells and b in sa_cells else "CR_only" if b in cr_cells else "SA_only" for b in union],
    index=union,
)
summary = group.value_counts().reindex(GROUP_ORDER, fill_value=0).rename("n_barcodes").to_frame()
summary.loc["CR_total"] = len(cr_cells)
summary.loc["SA_total"] = len(sa_cells)
summary.loc["jaccard"] = round(len(cr_cells & sa_cells) / len(union), 4)
summary.to_csv(TABLES / "barcode_sets_summary.tsv", sep="\t")
print(summary)

# --- 2. raw universes ----------------------------------------------------------
cr_umi, sa_umi = umi_per_bc(cr), umi_per_bc(sa)
uni = pd.DataFrame(
    {
        "n_barcodes": [len(cr_umi), len(sa_umi), len(set(cr_umi.index) & set(sa_umi.index))],
        "n_umi_ge_10": [(cr_umi >= 10).sum(), (sa_umi >= 10).sum(), np.nan],
        "n_umi_ge_100": [(cr_umi >= 100).sum(), (sa_umi >= 100).sum(), np.nan],
        "n_umi_ge_500": [(cr_umi >= 500).sum(), (sa_umi >= 500).sum(), np.nan],
        "total_umi": [cr_umi.sum(), sa_umi.sum(), np.nan],
    },
    index=["cellranger_raw", "simpleaf_raw", "intersection"],
)
only_cr = cr_umi[~cr_umi.index.isin(sa_umi.index)]
uni.loc["cellranger_only_tail", ["n_barcodes", "total_umi"]] = [len(only_cr), only_cr.sum()]
uni["max_umi"] = [cr_umi.max(), sa_umi.max(), np.nan, only_cr.max()]
uni.to_csv(TABLES / "raw_universe_summary.tsv", sep="\t")
print(uni)


# --- 3. per-barcode metrics from both quantifications -------------------------
def metrics(a, bcs, prefix):
    sub = a[bcs]
    X = sub.X
    umi = np.asarray(X.sum(1)).ravel()
    out = pd.DataFrame(
        {
            f"{prefix}_umi": umi,
            f"{prefix}_genes": np.asarray((X > 0).sum(1)).ravel(),
            f"{prefix}_pct_mt": 100 * np.asarray(X[:, mito_mask(sub)].sum(1)).ravel() / np.maximum(umi, 1),
        },
        index=bcs,
    )
    return out


def rank(u):
    return u.rank(ascending=False, method="min").astype(int)


m = pd.concat([metrics(cr, union, "cr"), metrics(sa, union, "sa")], axis=1)
m.insert(0, "group", group)
s_sub = sa[union]
m["sa_unspliced_frac"] = np.asarray(s_sub.layers["unspliced"].sum(1)).ravel() / np.maximum(m["sa_umi"], 1)
m["sa_ambiguous_frac"] = np.asarray(s_sub.layers["ambiguous"].sum(1)).ravel() / np.maximum(m["sa_umi"], 1)
m["cr_rank"] = rank(cr_umi).reindex(union).values
m["sa_rank"] = rank(sa_umi).reindex(union).values
m["umi_log2_ratio_sa_vs_cr"] = np.log2((m["sa_umi"] + 1) / (m["cr_umi"] + 1))
q = load_sa_qcatch_obs()
for c in ("initial_filtered_cell", "potential_non_ambient_cell", "non_ambient_pvalue"):
    m[f"qcatch_{c}"] = q[c].reindex(union).values
m.index.name = "barcode"
m.sort_values(["group", "cr_umi"], ascending=[True, False]).to_csv(TABLES / "union_cells_metrics.tsv", sep="\t")

print("\nDisputed barcodes:")
with pd.option_context("display.width", 250, "display.max_columns", 30):
    print(m[m.group != "shared"].round(3))


# --- 4. figures -----------------------------------------------------------------
def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, which="major", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), gridspec_kw={"width_ratios": [1, 1, 1.1]})
for ax, (name, u) in zip(axes[:2], [("Cell Ranger", cr_umi), ("Simpleaf", sa_umi)]):
    srt = np.sort(u.values)[::-1]
    srt = srt[srt > 0]
    ax.plot(np.arange(1, len(srt) + 1), srt, color=MUTED, linewidth=2, zorder=1)
    r = rank(u)
    for g in GROUP_ORDER:
        bcs = m.index[m.group == g]
        if g == "shared" or len(bcs) == 0:
            continue
        ax.scatter(r.reindex(bcs), u.reindex(bcs), s=64, marker=MARKERS[g], color=COLORS[g],
                   edgecolor="#fcfcfb", linewidth=2, zorder=3, label=f"{LABELS[g]} (n={len(bcs)})")
    for y, lab in [(500, "500 UMI (min. for EmptyDrops)")]:
        ax.axhline(y, color=MUTED, linewidth=1, linestyle=(0, (3, 3)))
        ax.text(1.2, y * 1.15, lab, color=MUTED, fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title(f"{name} raw matrix ({len(u):,} barcodes)", fontsize=11, color=INK, loc="left")
    ax.set_xlabel("Barcode rank", color=INK)
    style(ax)
axes[0].set_ylabel("UMI per barcode", color=INK)
axes[1].legend(frameon=False, fontsize=9, loc="lower left")

# Zoom on the knee: both curves on one axis, disputed barcodes numbered (same order as the tables).
ax = axes[2]
for name, u, ls in [("Cell Ranger", cr_umi, "-"), ("Simpleaf", sa_umi, (0, (4, 2)))]:
    srt = np.sort(u.values)[::-1]
    ax.plot(np.arange(1, len(srt) + 1), srt, color=INK if name == "Cell Ranger" else MUTED, linewidth=2,
            linestyle=ls, label=name, zorder=1)
dis = m[m.group != "shared"].sort_index()
for k, (b, row) in enumerate(dis.iterrows(), 1):
    ax.plot([row.cr_rank, row.sa_rank], [row.cr_umi, row.sa_umi], color=COLORS["SA_only"], linewidth=1, zorder=2)
    ax.scatter([row.cr_rank], [row.cr_umi], s=40, marker="o", facecolor="#fcfcfb", edgecolor=COLORS["SA_only"],
               linewidth=2, zorder=3)
    ax.scatter([row.sa_rank], [row.sa_umi], s=56, marker="D", color=COLORS["SA_only"], edgecolor="#fcfcfb",
               linewidth=1.5, zorder=4)
    ax.annotate(str(k), (row.sa_rank, row.sa_umi), xytext=(6, 2), textcoords="offset points", fontsize=8, color=INK)
ax.axhline(500, color=MUTED, linewidth=1, linestyle=(0, (3, 3)))
ax.text(1105, 520, "500 UMI", color=MUTED, fontsize=8)
ax.set_xlim(1100, 1320); ax.set_ylim(150, 3000); ax.set_yscale("log")
ax.set_title("Knee region: open circle = Cell Ranger, diamond = Simpleaf", fontsize=10, color=INK, loc="left")
ax.set_xlabel("Barcode rank", color=INK)
ax.legend(frameon=False, fontsize=9, loc="upper right")
style(ax)
fig.tight_layout()
fig.savefig(FIGURES / "01_knee_plots.png", dpi=200, facecolor="#fcfcfb")

fig, ax = plt.subplots(figsize=(5.6, 5.2))
for g in GROUP_ORDER:
    d = m[m.group == g]
    if len(d) == 0:
        continue
    ax.scatter(d.cr_umi, d.sa_umi, s=10 if g == "shared" else 64, marker=MARKERS[g], color=COLORS[g],
               alpha=0.5 if g == "shared" else 1, edgecolor="none" if g == "shared" else "#fcfcfb",
               linewidth=0 if g == "shared" else 2, label=f"{LABELS[g]} (n={len(d)})", zorder=2 if g == "shared" else 3)
lim = [300, m[["cr_umi", "sa_umi"]].max().max() * 1.3]
ax.plot(lim, lim, color=MUTED, linewidth=1, linestyle=(0, (3, 3)), zorder=1)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("UMI per barcode, Cell Ranger", color=INK)
ax.set_ylabel("UMI per barcode, Simpleaf (S+U+A)", color=INK)
ax.set_title("Same barcode, two quantifications", fontsize=11, color=INK, loc="left")
style(ax)
ax.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig(FIGURES / "01_umi_concordance.png", dpi=200, facecolor="#fcfcfb")

sh = m[m.group == "shared"]
print(f"\nShared cells: median log2(SA/CR) UMI = {sh.umi_log2_ratio_sa_vs_cr.median():.3f}; "
      f"Spearman r = {sh.cr_umi.corr(sh.sa_umi, method='spearman'):.4f}")
