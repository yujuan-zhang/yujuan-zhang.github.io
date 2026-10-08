"""Step 4: how much do cell-type proportions move for reasons unrelated to cell calling?

The same 1,228 barcodes are clustered and annotated with the step-3 workflow, varying only
  - the count matrix (Cell Ranger vs Simpleaf) and
  - the random seed of PCA / neighbours / Leiden / UMAP (0-4).
This gives the noise floor against which the cell-calling effect (<= 0.17 percentage points) is judged.

Output: tables/annotation_stability.tsv (one row per counts x seed x cell type)
"""

import numpy as np
import pandas as pd
import scanpy as sc

from common import TABLES, load_all

sc.settings.verbosity = 0
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
MIN_MARKER_SCORE, LOW_QUALITY = 0.3, "Low quality (no markers)"

cr, sa, cr_cells, sa_cells = load_all()
union = sorted(cr_cells | sa_cells)
rows = []
for label, raw in [("cellranger_counts", cr), ("simpleaf_counts", sa)]:
    base = raw[union].copy()
    base.var_names = base.var["gene_symbol"].astype(str).values
    base.var_names_make_unique()
    sc.pp.normalize_total(base, target_sum=1e4)
    sc.pp.log1p(base)
    sc.pp.highly_variable_genes(base, n_top_genes=2000, flavor="seurat")
    z = sc.pp.scale(base, max_value=10, copy=True)
    for seed in range(5):
        a = base.copy()
        sc.tl.pca(a, n_comps=30, mask_var="highly_variable", random_state=seed)
        sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=seed)
        sc.tl.leiden(a, resolution=0.8, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        lab = {}
        for cl in a.obs.leiden.unique():
            cells = (a.obs.leiden == cl).values
            scores = {ct: np.mean([z[cells, g].X.mean() for g in gs if g in z.var_names]) for ct, gs in MARKERS.items()}
            best = max(scores, key=scores.get)
            lab[cl] = best if scores[best] >= MIN_MARKER_SCORE else LOW_QUALITY
        types = a.obs.leiden.map(lab).astype(str)
        for setname, cells in [("Cell Ranger cells", cr_cells), ("Simpleaf+QCatch cells", sa_cells)]:
            vc = types[a.obs_names.isin(list(cells))].value_counts(normalize=True) * 100
            for ct, pct in vc.items():
                rows.append({"counts": label, "seed": seed, "cell_set": setname, "cell_type": ct, "pct": pct})
        print(label, seed, types.value_counts().to_dict(), flush=True)

st = pd.DataFrame(rows)
st.to_csv(TABLES / "annotation_stability.tsv", sep="\t", index=False)
full = st.pivot_table(index=["counts", "seed", "cell_set"], columns="cell_type", values="pct", fill_value=0)

# Three sources of variation, each as the max |change| in percentage points over cell types:
calling = (full.xs("Simpleaf+QCatch cells", level="cell_set") - full.xs("Cell Ranger cells", level="cell_set")).abs()
cr_set = full.xs("Cell Ranger cells", level="cell_set")
seed_range = cr_set.groupby(level="counts").agg(lambda s: s.max() - s.min())
counts_effect = (cr_set.xs("simpleaf_counts", level="counts") - cr_set.xs("cellranger_counts", level="counts")).abs()
print("\nCell calling (same counts, same seed): max |diff| per run =", calling.max(axis=1).round(2).tolist())
print("Seed only (same counts, same cell set): range per cell type")
print(seed_range.round(2).T)
print("Counts only (same seed, same cell set): max |diff| per seed =", counts_effect.max(axis=1).round(2).tolist())
summary = pd.DataFrame({
    "source_of_variation": ["cell calling (CR vs SA cell set)", "random seed (0-4)", "count matrix (CR vs SA)"],
    "max_abs_change_pp": [calling.values.max(), seed_range.values.max(), counts_effect.values.max()],
    "median_over_runs_of_max_change_pp": [np.median(calling.max(axis=1)), np.median(seed_range.max(axis=1)),
                                          np.median(counts_effect.max(axis=1))],
}).round(3)
summary.to_csv(TABLES / "annotation_stability_summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))
