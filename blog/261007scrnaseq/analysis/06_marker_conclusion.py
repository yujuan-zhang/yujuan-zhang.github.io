"""Step 6: does the cell-set difference change a concrete marker conclusion?

Cell types are those assigned in step 3 (Cell Ranger counts, identical for both cell sets).
For each cell type that gains disputed cells (cDC, Platelet, FCGR3A Mono), marker genes are
called with the same Wilcoxon test (type vs all other cells) on:
  - the Cell Ranger cell set (1,221) and
  - the Simpleaf + QCatch cell set (1,228),
both with and without the standard QC filter (>= 200 genes, <= 20 % mitochondrial).

Outputs: tables/marker_conclusion_top.tsv (top-10 genes per run), tables/marker_conclusion_genes.tsv
(canonical markers: log2FC, % expressing, adjusted p in each run).
"""

import numpy as np
import pandas as pd
import scanpy as sc

from common import TABLES, load_all

sc.settings.verbosity = 0
TYPES = {"cDC": ["FCER1A", "CD1C", "CLEC10A", "CST3"],
         "Platelet": ["PPBP", "PF4", "GNG11"],
         "FCGR3A Mono": ["FCGR3A", "MS4A7", "LST1"]}
TOP = 10

cr, sa, cr_cells, sa_cells = load_all()
ann = pd.read_csv(TABLES / "union_cells_annotated.tsv", sep="\t", index_col=0)
base = cr[list(ann.index)].copy()
base.var_names = base.var["gene_symbol"].astype(str).values
base.var_names_make_unique()
base.obs["cell_type"] = ann["cr_cell_type"].values
base.obs["qc"] = ann["cr_passes_std_qc"].values
sc.pp.normalize_total(base, target_sum=1e4)
sc.pp.log1p(base)

tops, genes = [], []
for setname, cells in [("Cell Ranger", cr_cells), ("Simpleaf+QCatch", sa_cells)]:
    for qc in (False, True):
        a = base[base.obs_names.isin(list(cells)) & (base.obs.qc.values if qc else True)].copy()
        for ct, markers in TYPES.items():
            a.obs["is_type"] = np.where(a.obs.cell_type == ct, ct, "rest")
            n = int((a.obs.is_type == ct).sum())
            sc.tl.rank_genes_groups(a, "is_type", groups=[ct], reference="rest", method="wilcoxon", pts=True)
            df = sc.get.rank_genes_groups_df(a, group=ct)
            df = df.sort_values("scores", ascending=False).reset_index(drop=True)
            tops.append({"cell_set": setname, "after_qc": qc, "cell_type": ct, "n_cells": n,
                         "top_genes": ", ".join(df.names[:TOP])})
            for g in markers:
                r = df[df.names == g].iloc[0]
                genes.append({"cell_set": setname, "after_qc": qc, "cell_type": ct, "n_cells": n, "gene": g,
                              "rank": int(df.index[df.names == g][0]) + 1, "log2FC": round(r.logfoldchanges, 2),
                              "pct_in_type": round(100 * r.pct_nz_group, 1),
                              "pct_in_rest": round(100 * r.pct_nz_reference, 1), "padj": r.pvals_adj})

top = pd.DataFrame(tops)
gen = pd.DataFrame(genes)
top.to_csv(TABLES / "marker_conclusion_top.tsv", sep="\t", index=False)
gen.to_csv(TABLES / "marker_conclusion_genes.tsv", sep="\t", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 140)
for (ct, qc), d in top.groupby(["cell_type", "after_qc"]):
    a, b = d.iloc[0], d.iloc[1]
    sa_set, cr_set = set(b.top_genes.split(", ")), set(a.top_genes.split(", "))
    print(f"\n{ct}, after_qc={qc}: n {a.n_cells} vs {b.n_cells}; top-{TOP} overlap {len(sa_set & cr_set)}/{TOP}")
    print("  Cell Ranger    :", a.top_genes)
    print("  Simpleaf+QCatch:", b.top_genes)
w = gen.pivot_table(index=["cell_type", "gene", "after_qc"], columns="cell_set",
                    values=["rank", "log2FC", "pct_in_type", "padj"], aggfunc="first")
print("\n", w.to_string(float_format=lambda x: f"{x:.3g}"))
