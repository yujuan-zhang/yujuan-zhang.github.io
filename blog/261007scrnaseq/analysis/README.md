# Comparing the cell calls of Cell Ranger and Simpleaf

The input is the output of the AWS run `../scrna-results/20261005T210235Z/` (nf-core/scrnaseq 4.2.0, PBMC 1k v3).
Paths are set in [common.py](common.py).

## Environment

```bash
/path/to/python3.11 -m venv ../.venv
../.venv/bin/pip install "numba==0.61.2" "llvmlite==0.44.0"   # newer releases have no Intel-Mac wheels
../.venv/bin/pip install -r requirements.txt
```

QCatch is pinned to **0.2.12**, the version that ran in the pipeline (recorded in `uns/qc_info/version`), so step 2 runs the same cell-calling code.

## Run order

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `01_barcode_sets.py` | Three barcode groups; per-barcode metrics computed from **both** raw matrices; raw barcode universes | `tables/barcode_sets_summary.tsv`, `tables/union_cells_metrics.tsv`, `tables/raw_universe_summary.tsv`, `figures/01_*.png` |
| `02_cell_calling_decomposition.py` | Re-runs QCatch's cell calling with swapped inputs (counts × barcode universe × simulation count × random seed) | `tables/cell_calling_conditions_*.tsv`, `tables/cell_calling_barcodes_*.tsv` |
| `03_disputed_cells.py` | Cell type, markers, doublet score, depth-matched comparison with ambient RNA, cell-type proportions | `tables/union_cells_annotated.tsv`, `tables/disputed_depth_matched.tsv`, `tables/celltype_proportions.tsv`, `figures/03_*.png` |
| `04_annotation_stability.py` | Noise floor for proportions: random seed and count matrix vs cell calling | `tables/annotation_stability*.tsv` |
| `05_dhfr_locus.py` | Read-level check of DHFR and the chr8 lncRNA in the Cell Ranger BAM (read remotely from S3) | `tables/dhfr_locus_reads.tsv` |
| `06_marker_conclusion.py` | Marker genes for cDC, platelets and FCGR3A monocytes under each cell set (same labels, same Wilcoxon test) | `tables/marker_conclusion_top.tsv`, `tables/marker_conclusion_genes.tsv` |
| `07_impact_figures.py` | Figures for the effect on conclusions: per-type proportion change by source; cDC markers under both cell sets | `tables/proportion_impact_by_type.tsv`, `figures/07_*.png` |

Step 2 takes about 4–12 minutes per condition at 10,000 simulations and roughly ten times longer at 100,000. Each condition can run as a separate process:

```bash
python 02_cell_calling_decomposition.py --conditions A --sims 10000 --seed 42
python 02_cell_calling_decomposition.py --conditions B,C --sims 10000 --seed 42
python 02_cell_calling_decomposition.py --conditions D,E --sims 10000 --seed 42
python 02_cell_calling_decomposition.py --conditions A,B --sims 100000 --seed 42
for s in 1 2 3; do python 02_cell_calling_decomposition.py --conditions A,B --sims 10000 --seed $s; done
```

| Condition | Counts | Barcode universe |
| --- | --- | --- |
| A | Simpleaf (S+U+A) | Simpleaf raw, 72,612 |
| B | Cell Ranger | Cell Ranger raw, 329,735 |
| C | Cell Ranger | Barcodes shared by both, 71,929 |
| D | Simpleaf | Barcodes shared by both, 71,929 |
| E | Simpleaf | Simpleaf raw plus the Cell Ranger-only tail, 330,418 |
| F | Simpleaf, without DHFR and ENSG00000247134 | Simpleaf raw |
| G | Cell Ranger, without DHFR and ENSG00000247134 | Cell Ranger raw |

`05_dhfr_locus.py` reads only the DHFR and chr8 lncRNA regions from the Cell Ranger BAM in S3, through presigned URLs (usage is in the script header). It needs `pysam`.

The findings are in [FINDINGS.md](FINDINGS.md). `cache/` and `logs/` are not tracked.
