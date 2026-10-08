"""Shared paths and loaders for the Cell Ranger vs Simpleaf comparison.

All matrices are returned as AnnData with barcodes (no "-1" suffix) as obs_names
and Ensembl gene IDs as var_names, in the same gene order for both routes.
Counts are integers. For Simpleaf, X = spliced + unspliced + ambiguous, which is
the intron-inclusive count that matches Cell Ranger's default.
"""

from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "scrna-results" / "20261005T210235Z" / "results"
CR_OUTS = RUN / "cellranger" / "cellranger" / "count" / "PBMC1K" / "outs"
SA_DIR = RUN / "simpleaf" / "simpleaf" / "PBMC1K"
SA_RAW = SA_DIR / "simpleaf_quant" / "af_quant" / "alevin" / "quants.h5ad"
SA_QCATCH = SA_DIR / "qcatch" / "PBMC1K_filtered_quants.h5ad"

OUT = ROOT / "analysis"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
CACHE = OUT / "cache"
for d in (TABLES, FIGURES, CACHE):
    d.mkdir(exist_ok=True)

GROUP_ORDER = ["shared", "CR_only", "SA_only"]


def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def read_10x_h5(path: Path) -> ad.AnnData:
    """Read a Cell Ranger feature-barcode h5 as barcodes x genes (integer counts)."""
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        bcs = np.char.partition(_dec(m["barcodes"][:]).astype(str), "-")[:, 0]
        genes = _dec(m["features/id"][:])
        symbols = _dec(m["features/name"][:])
        shape = tuple(m["shape"][:])  # (genes, barcodes)
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape)
    a = ad.AnnData(X.T.tocsr().astype(np.int64))
    a.obs_names = bcs
    a.var_names = genes
    a.var["gene_symbol"] = symbols
    return a


def load_cr_raw() -> ad.AnnData:
    return read_10x_h5(CR_OUTS / "raw_feature_bc_matrix.h5")


def load_cr_cells() -> set[str]:
    with h5py.File(CR_OUTS / "filtered_feature_bc_matrix.h5", "r") as f:
        return set(np.char.partition(_dec(f["matrix/barcodes"][:]).astype(str), "-")[:, 0])


def load_sa_raw(gene_order=None) -> ad.AnnData:
    """Simpleaf unfiltered-permit-list matrix; X = S+U+A (checked equal to the layer sum)."""
    a = ad.read_h5ad(SA_RAW)
    total = sum(a.layers[k] for k in ("spliced", "unspliced", "ambiguous"))
    assert abs(total - a.X).max() == 0, "Simpleaf X is not S+U+A"
    a.X = sp.csr_matrix(a.X).astype(np.int64)
    for k in ("spliced", "unspliced", "ambiguous"):
        a.layers[k] = sp.csr_matrix(a.layers[k]).astype(np.int64)
    a.obs_names = a.obs["barcodes"].astype(str).values
    a.var_names = a.var["gene_id"].astype(str).values
    if gene_order is not None:
        assert set(gene_order) == set(a.var_names), "gene sets differ between routes"
        a = a[:, list(gene_order)].copy()
    return a


def load_sa_qcatch_obs() -> pd.DataFrame:
    """Per-barcode QCatch cell-calling fields for the 1,228 retained cells."""
    a = ad.read_h5ad(SA_QCATCH, backed="r")
    obs = a.obs.copy()
    obs.index = obs["barcodes"].astype(str).values
    return obs


def load_all():
    """Return (cr_raw, sa_raw, cr_cells, sa_cells) with identical gene order."""
    cr = load_cr_raw()
    sa = load_sa_raw(gene_order=cr.var_names)
    sa.var["gene_symbol"] = cr.var["gene_symbol"].values
    cr_cells = load_cr_cells()
    sa_cells = set(load_sa_qcatch_obs().index)
    return cr, sa, cr_cells, sa_cells


def umi_per_bc(a: ad.AnnData) -> pd.Series:
    return pd.Series(np.asarray(a.X.sum(axis=1)).ravel(), index=a.obs_names)


def mito_mask(a: ad.AnnData) -> np.ndarray:
    return a.var["gene_symbol"].astype(str).str.upper().str.startswith("MT-").values
