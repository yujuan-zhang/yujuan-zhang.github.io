"""Step 3: what are the disputed barcodes, and do they change cell-type proportions?

All 1,228 called barcodes (union) are analysed together, once with Cell Ranger counts and
once with Simpleaf counts, using the same scanpy workflow and the same marker-based annotation.

Outputs
  tables/cluster_annotation_<counts>.tsv  cluster -> cell type, with marker evidence
  tables/union_cells_annotated.tsv        per barcode: cell type (both quantifications), scrublet,
                                          ambient vs cell-type similarity, standard QC pass/fail
  tables/celltype_proportions.tsv         counts and % per cell type under each tool's cell set
  figures/03_umap_disputed.png, 03_qc_disputed.png, 03_marker_dotplot.png
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from common import FIGURES, TABLES, load_all, mito_mask

sc.settings.verbosity = 0
COLORS = {"shared": "#2a78d6", "SA_only": "#1baf7a", "CR_only": "#eb6834"}
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SEED = 0

# Canonical PBMC markers. A cluster takes the label whose markers have the highest mean scaled expression.
MARKERS = {
    "CD4 T": ["CD3E", "IL7R", "CCR7", "LEF1", "CD4"],
    "CD8 T": ["CD3E", "CD8A", "CD8B"],
    "NK": ["NKG7", "GNLY", "KLRD1", "PRF1"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "CD14 Mono": ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A Mono": ["FCGR3A", "MS4A7", "LST1"],
    "cDC": ["FCER1A", "CST3", "CD1C"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C"],
    "Platelet": ["PPBP", "PF4"],
}
MIN_MARKER_SCORE = 0.3  # mean scaled marker expression a cluster needs to receive a cell-type label
LOW_QUALITY = "Low quality (no markers)"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)

cr, sa, cr_cells, sa_cells = load_all()
union = sorted(cr_cells | sa_cells)
group = pd.Series(["shared" if b in cr_cells and b in sa_cells else "CR_only" if b in cr_cells else "SA_only"
                   for b in union], index=union)
metrics = pd.read_csv(TABLES / "union_cells_metrics.tsv", sep="\t", index_col=0)

# Ambient profile: barcodes in the EmptyDrops background range of the Cell Ranger raw matrix (ranks 45k-90k).
cr_umi = np.asarray(cr.X.sum(1)).ravel()
amb_idx = np.argsort(cr_umi)[::-1][45000:90000]
ambient = np.asarray(cr.X[amb_idx].sum(0)).ravel().astype(float)
ambient /= ambient.sum()


def analyse(raw, label):
    a = raw[union].copy()
    a.var_names = a.var["gene_symbol"].astype(str).values
    a.var_names_make_unique()
    a.obs["group"] = group.values
    a.layers["counts"] = a.X.copy()
    a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    sc.pp.scrublet(a, random_state=SEED)  # on raw counts, before normalisation
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
    sc.tl.pca(a, n_comps=30, mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.leiden(a, resolution=0.8, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(a, random_state=SEED)

    # Marker-based cluster annotation (scaled expression, cluster means).
    z = sc.pp.scale(a, max_value=10, copy=True)
    rows = []
    for cl in sorted(a.obs.leiden.unique(), key=int):
        cells = (a.obs.leiden == cl).values
        scores = {ct: float(np.mean([z[cells, g].X.mean() for g in genes if g in z.var_names]))
                  for ct, genes in MARKERS.items()}
        best = max(scores, key=scores.get)
        second = sorted(scores.values())[-2]
        rows.append({"cluster": cl, "n_cells": int(cells.sum()), "cell_type": best,
                     "best_score": round(scores[best], 2), "margin_over_2nd": round(scores[best] - second, 2),
                     **{f"score_{k}": round(v, 2) for k, v in scores.items()}})
    ann = pd.DataFrame(rows)
    # A cluster with no positive marker signal is not a cell type: in this sample it is the
    # high-mitochondrial, few-gene cluster that would otherwise win "Platelet" by a 0.07 margin.
    ann.loc[ann.best_score < MIN_MARKER_SCORE, "cell_type"] = LOW_QUALITY
    ann.to_csv(TABLES / f"cluster_annotation_{label}.tsv", sep="\t", index=False)
    a.obs["cell_type"] = a.obs.leiden.map(dict(zip(ann.cluster, ann.cell_type))).astype(str)

    # Similarity of each cell to the ambient profile vs to its own cell type's centroid (log-CPM, Pearson).
    logcpm = a.X.toarray()
    amb_log = np.log1p(ambient * 1e4)
    cent = {ct: np.asarray(logcpm[(a.obs.cell_type == ct).values & (a.obs.group == "shared").values].mean(0)).ravel()
            for ct in a.obs.cell_type.unique()}
    a.obs["r_ambient"] = [np.corrcoef(logcpm[i], amb_log)[0, 1] for i in range(a.n_obs)]
    a.obs["r_own_type"] = [np.corrcoef(logcpm[i], cent[ct])[0, 1] for i, ct in enumerate(a.obs.cell_type)]

    # kNN label from shared cells only (does the cell sit inside a population of agreed cells?)
    conn = a.obsp["connectivities"].tocsr()
    shared_mask = (a.obs.group == "shared").values
    knn_lab, knn_frac = [], []
    for i in range(a.n_obs):
        nb = conn[i].indices[shared_mask[conn[i].indices]]
        labs = a.obs.cell_type.values[nb]
        if len(labs) == 0:
            knn_lab.append("none"); knn_frac.append(np.nan); continue
        v, c = np.unique(labs, return_counts=True)
        knn_lab.append(v[c.argmax()]); knn_frac.append(c.max() / c.sum())
    a.obs["knn_type_from_shared"] = knn_lab
    a.obs["knn_agreement"] = knn_frac
    return a


res = {lab: analyse(raw, lab) for lab, raw in [("cellranger_counts", cr), ("simpleaf_counts", sa)]}
A, B = res["cellranger_counts"], res["simpleaf_counts"]

# Standard downstream QC that many tutorials apply (fixed thresholds, stated up front).
QC = {"min_genes": 200, "max_pct_mt": 20}
out = pd.DataFrame(index=union)
out["group"] = group
for lab, a in res.items():
    tag = "cr" if lab.startswith("cell") else "sa"
    o = a.obs
    out[f"{tag}_cell_type"] = o.cell_type.values
    out[f"{tag}_knn_type_from_shared"] = o.knn_type_from_shared.values
    out[f"{tag}_knn_agreement"] = o.knn_agreement.round(3).values
    out[f"{tag}_r_ambient"] = o.r_ambient.round(3).values
    out[f"{tag}_r_own_type"] = o.r_own_type.round(3).values
    out[f"{tag}_doublet_score"] = o.doublet_score.round(3).values
    out[f"{tag}_predicted_doublet"] = o.predicted_doublet.values
    out[f"{tag}_passes_std_qc"] = ((o.n_genes_by_counts >= QC["min_genes"]) & (o.pct_counts_mt <= QC["max_pct_mt"])).values
out = out.join(metrics.drop(columns="group"))
out.index.name = "barcode"
out.to_csv(TABLES / "union_cells_annotated.tsv", sep="\t")

agree = (out.cr_cell_type == out.sa_cell_type).mean()
print(f"Cell-type label agreement between the two quantifications (all 1,228): {agree:.2%}")
with pd.option_context("display.width", 250, "display.max_columns", 40):
    cols = ["cr_cell_type", "sa_cell_type", "cr_knn_type_from_shared", "cr_knn_agreement", "cr_r_ambient",
            "cr_r_own_type", "cr_doublet_score", "cr_passes_std_qc", "sa_passes_std_qc", "cr_umi", "cr_genes", "cr_pct_mt"]
    print(out.loc[out.group != "shared", cols])
    sh = out[out.group == "shared"]
    print("\nShared cells, reference distribution (5th / 50th / 95th percentile):")
    print(sh[["cr_r_ambient", "cr_r_own_type", "cr_doublet_score", "cr_umi", "cr_genes", "cr_pct_mt"]]
          .quantile([0.05, 0.5, 0.95]).round(3))

# --- depth-matched check: does each disputed barcode look like a cell or like ambient RNA? -------------
# Correlations depend on depth, so each disputed barcode (n UMIs) is compared with
#   (i) agreed cells of its own type, downsampled to n UMIs, and
#  (ii) n-UMI multinomial draws from the ambient profile (what an empty droplet would look like).
N_REF = 200
rng = np.random.default_rng(SEED)
counts = A.layers["counts"].tocsr()
amb_log = np.log1p(ambient * 1e4)
# A.X holds log1p(CP10k), the same transform two_r applies, so centroids and test vectors match.
cent ={ct: np.asarray(A.X[(A.obs.cell_type == ct).values & (A.obs.group == "shared").values].mean(0)).ravel()
        for ct in A.obs.cell_type.unique()}


def two_r(vec, centroid):
    v = np.log1p(vec / max(vec.sum(), 1) * 1e4)
    return np.corrcoef(v, centroid)[0, 1], np.corrcoef(v, amb_log)[0, 1]


dm_rows, dm_points = [], []
for b in out.index[out.group != "shared"]:
    i = A.obs_names.get_loc(b)
    n = int(counts[i].sum())
    ct = A.obs.cell_type.iloc[i]
    obs_type, obs_amb = two_r(counts[i].toarray().ravel(), cent[ct])
    pool = np.flatnonzero((A.obs.group == "shared").values & (A.obs.cell_type == ct).values
                          & (np.asarray(counts.sum(1)).ravel() >= n))
    pick = rng.choice(pool, size=min(N_REF, len(pool)), replace=False)
    ds = np.array([two_r(rng.multivariate_hypergeometric(counts[j].toarray().ravel(), n), cent[ct]) for j in pick])
    em = np.array([two_r(rng.multinomial(n, ambient), cent[ct]) for _ in range(N_REF)])
    # Distance (in r_own_type units) to the median of each reference cloud.
    d_cell, d_amb = abs(obs_type - np.median(ds[:, 0])), abs(obs_type - np.median(em[:, 0]))
    dm_rows.append({"barcode": b, "cell_type": ct, "umi": n, "n_reference_cells": len(pick),
                    "r_own_type": round(obs_type, 3), "r_ambient": round(obs_amb, 3),
                    "ref_cells_r_own_type_median": round(np.median(ds[:, 0]), 3),
                    "ref_cells_r_own_type_5pct": round(np.percentile(ds[:, 0], 5), 3),
                    "ambient_draws_r_own_type_median": round(np.median(em[:, 0]), 3),
                    "ambient_draws_r_own_type_95pct": round(np.percentile(em[:, 0], 95), 3),
                    "closer_to": "cells" if d_cell < d_amb else "ambient"})
    dm_points.append((b, ct, n, obs_type, obs_amb, ds, em))
dm = pd.DataFrame(dm_rows).set_index("barcode")
dm.to_csv(TABLES / "disputed_depth_matched.tsv", sep="\t")
print("\nDepth-matched comparison (Cell Ranger counts):")
with pd.option_context("display.width", 250, "display.max_columns", 20):
    print(dm)

fig, axes = plt.subplots(1, len(dm_points), figsize=(2.6 * len(dm_points), 3.0), sharex=True, sharey=True)
for ax, (b, ct, n, ot, oa, ds, em) in zip(np.atleast_1d(axes), dm_points):
    ax.scatter(em[:, 1], em[:, 0], s=6, color="#a3a29b", edgecolor="none", label="Ambient draws")
    ax.scatter(ds[:, 1], ds[:, 0], s=6, color=COLORS["shared"], alpha=0.6, edgecolor="none", label="Agreed cells, downsampled")
    ax.scatter([oa], [ot], s=70, marker="D", color=COLORS["SA_only"], edgecolor=SURFACE, linewidth=2, zorder=3,
               label="Disputed barcode")
    ax.set_title(f"{b[:8]}\n{ct}, {n} UMI", fontsize=8, color=INK, loc="left")
    ax.grid(True, color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    style(ax)
np.atleast_1d(axes)[0].set_ylabel("r with own cell type", color=INK, fontsize=9)
for ax in np.atleast_1d(axes):
    ax.set_xlabel("r with ambient", color=INK, fontsize=9)
np.atleast_1d(axes)[-1].legend(frameon=False, fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(FIGURES / "03_depth_matched.png", dpi=200, facecolor=SURFACE)

# --- proportions under each tool's cell set (labels from the same annotation run) --------------------
rows = []
for lab in ("cr", "sa"):
    types = out[f"{lab}_cell_type"]
    for setname, cells in [("Cell Ranger cells", cr_cells), ("Simpleaf+QCatch cells", sa_cells)]:
        for qc in (False, True):
            sel = out.index.isin(list(cells))
            if qc:
                sel &= out[f"{lab}_passes_std_qc"].values
            vc = types[sel].value_counts()
            for ct, n in vc.items():
                rows.append({"labels_from": f"{lab}_counts", "cell_set": setname, "after_std_qc": qc,
                             "cell_type": ct, "n": int(n), "pct": round(100 * n / sel.sum(), 3), "total": int(sel.sum())})
prop = pd.DataFrame(rows)
prop.to_csv(TABLES / "celltype_proportions.tsv", sep="\t", index=False)
wide = prop[prop.labels_from == "cr_counts"].pivot_table(index="cell_type", columns=["after_std_qc", "cell_set"], values="pct")
print("\nCell-type % (labels from Cell Ranger counts):")
print(wide.round(2))
for qc in (False, True):
    d = (wide[(qc, "Simpleaf+QCatch cells")] - wide[(qc, "Cell Ranger cells")]).abs().max()
    print(f"max |difference| in percentage points, after_std_qc={qc}: {d:.3f}")

# --- figures --------------------------------------------------------------------------------------------
fig, (ax, side) = plt.subplots(1, 2, figsize=(9.6, 6), gridspec_kw={"width_ratios": [3.2, 1]})
um = A.obsm["X_umap"]
types = A.obs.cell_type.values
ax.scatter(um[:, 0], um[:, 1], s=6, color="#c9c8c2", edgecolor="none", zorder=1)
for ct in sorted(set(types)):
    pts = um[types == ct]
    medoid = pts[np.argmin(((pts - np.median(pts, axis=0)) ** 2).sum(1))]  # a real point inside the cluster
    print(f"label {ct}: n={len(pts)} at {medoid.round(2)}")
    ax.annotate(ct, medoid, xytext=(0, -14), textcoords="offset points", fontsize=9, color=INK, ha="center",
                va="top", zorder=4, bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none", alpha=0.85))
dis_idx = [A.obs_names.get_loc(b) for b in dm.index]
OFFSETS = [(8, 6), (8, -12), (-14, 8), (8, 10), (-14, -12), (8, -12), (8, 8)]
ax.scatter(um[dis_idx, 0], um[dis_idx, 1], s=70, marker="D", color=COLORS["SA_only"], edgecolor=SURFACE,
           linewidth=2, zorder=5, label=f"Simpleaf + QCatch only (n={len(dis_idx)})")
for k, i in enumerate(dis_idx, 1):
    ax.annotate(str(k), um[i], xytext=OFFSETS[(k - 1) % len(OFFSETS)], textcoords="offset points", fontsize=9, color=INK,
                fontweight="bold", zorder=6)
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("1,228 called barcodes (UMAP, Cell Ranger counts)", fontsize=11, color=INK, loc="left")
for s_ in ax.spines.values():
    s_.set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper left")
side.axis("off")
side.text(0, 1, "Disputed barcodes", fontsize=10, color=INK, fontweight="bold", va="top")
for k, b in enumerate(dm.index, 1):
    side.text(0, 1 - 0.07 * k, f"{k}  {b[:8]}  {dm.loc[b, 'cell_type']}, {dm.loc[b, 'umi']} UMI",
              fontsize=8, color=INK, va="top")
fig.tight_layout(); fig.savefig(FIGURES / "03_umap_disputed.png", dpi=200, facecolor=SURFACE)

qc_cols = [("cr_umi", "UMI (Cell Ranger)", True), ("cr_genes", "Genes detected", True),
           ("cr_pct_mt", "% mitochondrial", False), ("sa_unspliced_frac", "Unspliced fraction (Simpleaf)", False),
           ("umi_log2_ratio_sa_vs_cr", "log2 UMI ratio, Simpleaf / Cell Ranger", False), ("cr_doublet_score", "Scrublet score", False)]
fig, axes = plt.subplots(1, len(qc_cols), figsize=(15, 3.6))
rng = np.random.default_rng(SEED)
for ax, (col, title, logy) in zip(axes, qc_cols):
    sh = out.loc[out.group == "shared", col]
    ax.scatter(rng.normal(0, 0.06, len(sh)), sh, s=4, color=COLORS["shared"], alpha=0.35, edgecolor="none")
    d = out.loc[out.group != "shared", col]
    ax.scatter(0.45 + np.linspace(-0.12, 0.12, len(d)), d, s=48, marker="D", color=COLORS["SA_only"], edgecolor=SURFACE, linewidth=1.5)
    ax.set_xticks([0, 0.45]); ax.set_xticklabels(["Both", "SA only"], fontsize=8)
    ax.set_xlim(-0.35, 0.75)
    if logy:
        ax.set_yscale("log")
    ax.set_title(title, fontsize=9, color=INK, loc="left")
    ax.grid(True, axis="y", color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    style(ax)
fig.suptitle("Disputed barcodes against the 1,221 agreed cells", x=0.01, ha="left", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(FIGURES / "03_qc_disputed.png", dpi=200, facecolor=SURFACE)

genes = [g for gs in MARKERS.values() for g in gs if g in A.var_names]
genes = list(dict.fromkeys(genes))
A.obs["type_or_disputed"] = np.where(A.obs.group == "shared", A.obs.cell_type, "SA_only: " + A.obs_names.str[:6])
dp = sc.pl.dotplot(A, genes, groupby="type_or_disputed", standard_scale="var", show=False, return_fig=True)
dp.savefig(FIGURES / "03_marker_dotplot.png", dpi=200, bbox_inches="tight", facecolor=SURFACE)
print("figures written")
