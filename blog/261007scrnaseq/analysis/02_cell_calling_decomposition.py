"""Step 2: why do the calls differ? Re-run one cell-calling implementation on swapped inputs.

Cell Ranger 10.0.0 (3' v3 GEX) and QCatch 0.2.12 implement the same two-step algorithm:
OrdMag initial calls, then EmptyDrops on non-cell barcodes with >= max(500, 1 + max ambient UMI)
UMIs, ambient = barcodes ranked 45,000-90,000 by UMI, BH FDR <= 0.001. Known implementation
differences: Monte Carlo simulations (CR 100,000; QCatch 10,000) and the random streams.

We run QCatch's own functions (installed qcatch==0.2.12) on:
  counts   : Cell Ranger vs Simpleaf (S+U+A) quantification
  universe : which barcodes are in the matrix (each tool's raw, or the shared intersection)
  num_sims : 10,000 (QCatch) vs 100,000 (Cell Ranger)

Usage: python 02_cell_calling_decomposition.py [--sims 10000] [--conditions A,B,...] [--seed 42]
Outputs: tables/cell_calling_conditions[_<sims>].tsv, tables/cell_calling_barcodes[_<sims>].tsv
"""

import argparse
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from qcatch.logger import setup_logger

logger = setup_logger("qcatch", False)
from qcatch.find_retained_cells import cell_calling as cc  # noqa: E402  (needs the logger first)
from qcatch.find_retained_cells.matrix import CountMatrix  # noqa: E402

from common import TABLES, load_all  # noqa: E402

CHEM = "10X_3p_v3"
OUTLIER_GENES = ["ENSG00000228716", "ENSG00000247134"]  # DHFR; lncRNA on chr8 (see FINDINGS.md)
_orig_sim = cc.simulate_multinomial_loglikelihoods


def call_cells(adata: ad.AnnData, num_sims: int, seed: int) -> dict:
    """QCatch internal_cell_calling, with the simulation count and RNG seed made explicit."""
    cc.RNG = np.random.default_rng(seed)  # QCatch seeds a module-level generator with 42
    cc.simulate_multinomial_loglikelihoods = lambda *a, **k: _orig_sim(*a, **{**k, "num_sims": num_sims})
    a = adata.copy()
    a.obs["barcodes"] = a.obs_names.values
    m = CountMatrix.from_anndata(a)
    filtered_bcs = cc.initial_filtering_OrdMag(m, CHEM, None)
    initial = [x.decode() if isinstance(x, bytes) else str(x) for x in filtered_bcs]
    # As in qcatch run_cell_calling: pass the un-decoded barcodes, which match matrix.bcs.
    res = cc.find_nonambient_barcodes(m, filtered_bcs, CHEM, None)
    assert res is not None, "EmptyDrops step returned None"

    # Recompute the quantities the algorithm uses internally, for reporting.
    umis = np.asarray(a.X.sum(1)).ravel()
    lo, hi = cc.compute_empty_drops_bounds(CHEM, None)
    order = np.argsort(umis)[::-1]
    empty = order[lo:hi]
    max_amb = int(umis[empty].max()) if len(empty) else 0
    return {
        "initial": set(initial),
        "eval": pd.DataFrame(
            {"pvalue": res.pvalues, "pvalue_adj": res.pvalues_adj, "nonambient": res.is_nonambient},
            index=res.eval_bcs,
        ),
        "n_barcodes": a.n_obs,
        "ambient_rank_range": f"{lo}-{min(hi, a.n_obs)}",
        "n_ambient_bcs": int(len(empty)),
        "max_ambient_umi": max_amb,
        "ambient_umi_range": f"{int(umis[empty].min()) if len(empty) else 0}-{max_amb}",
        "min_umi_to_test": max(cc.MIN_UMIS, 1 + max_amb),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sims", type=int, default=10000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--conditions", default="A,B,C,D,E")
    args = p.parse_args()

    cr, sa, cr_cells, sa_cells = load_all()
    shared_bcs = cr.obs_names.intersection(sa.obs_names)
    cr_tail = cr[~cr.obs_names.isin(sa.obs_names)]
    conditions = {
        "A": ("Simpleaf counts", "Simpleaf raw", lambda: sa),
        "B": ("Cell Ranger counts", "Cell Ranger raw", lambda: cr),
        "C": ("Cell Ranger counts", "shared barcodes", lambda: cr[shared_bcs]),
        "D": ("Simpleaf counts", "shared barcodes", lambda: sa[shared_bcs]),
        "E": ("Simpleaf counts", "Simpleaf raw + Cell Ranger-only tail",
              lambda: ad.AnnData(sp.vstack([sa.X, cr_tail.X]).tocsr(),
                                 obs=pd.DataFrame(index=list(sa.obs_names) + list(cr_tail.obs_names)),
                                 var=sa.var[[]])),
        # Causal test: drop the two loci where Simpleaf counts far more unspliced UMIs than Cell Ranger.
        "F": ("Simpleaf counts without DHFR and ENSG00000247134", "Simpleaf raw",
              lambda: sa[:, ~sa.var_names.isin(OUTLIER_GENES)]),
        "G": ("Cell Ranger counts without DHFR and ENSG00000247134", "Cell Ranger raw",
              lambda: cr[:, ~cr.var_names.isin(OUTLIER_GENES)]),
    }
    disputed = sorted(sa_cells - cr_cells) + sorted(cr_cells - sa_cells)
    rows, per_bc = [], {}
    for key in args.conditions.split(","):
        counts, universe, make = conditions[key]
        t0 = time.time()
        r = call_cells(make(), args.sims, args.seed)
        called = r["initial"] | set(r["eval"].index[r["eval"].nonambient])
        rows.append({
            "condition": key, "counts": counts, "universe": universe, "num_sims": args.sims, "seed": args.seed,
            "n_barcodes": r["n_barcodes"], "ambient_rank_range": r["ambient_rank_range"],
            "ambient_umi_range": r["ambient_umi_range"], "min_umi_to_test": r["min_umi_to_test"],
            "n_initial_ordmag": len(r["initial"]), "n_tested_emptydrops": len(r["eval"]),
            "n_added_emptydrops": int(r["eval"].nonambient.sum()), "n_called": len(called),
            "overlap_cellranger_official": len(called & cr_cells),
            "overlap_qcatch_official": len(called & sa_cells),
            "identical_to_cellranger": called == cr_cells, "identical_to_qcatch": called == sa_cells,
            "seconds": round(time.time() - t0, 1),
        })
        print(rows[-1], flush=True)
        status = {}
        for b in disputed:
            if b in r["initial"]:
                status[b] = "ordmag"
            elif b in r["eval"].index:
                e = r["eval"].loc[b]
                status[b] = f"{'CALLED' if e.nonambient else 'not called'} (p={e.pvalue:.2e}, padj={e.pvalue_adj:.2e})"
            else:
                status[b] = "not tested (UMI below threshold or absent)"
        per_bc[key] = status
        # Barcodes whose call differs from the official results, beyond the 7 known.
        extra = (called ^ cr_cells) - set(disputed)
        if extra:
            print(f"  condition {key}: {len(extra)} other barcodes differ from Cell Ranger: {sorted(extra)[:10]}")

    suffix = f"_{args.conditions.replace(',', '')}_sims{args.sims}_seed{args.seed}"
    pd.DataFrame(rows).to_csv(TABLES / f"cell_calling_conditions{suffix}.tsv", sep="\t", index=False)
    pd.DataFrame(per_bc).rename_axis("barcode").to_csv(TABLES / f"cell_calling_barcodes{suffix}.tsv", sep="\t")


if __name__ == "__main__":
    main()
