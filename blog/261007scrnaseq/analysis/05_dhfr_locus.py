"""Step 5: read-level check of the DHFR locus in the Cell Ranger BAM.

The BAM stays in S3; only the needed region is read through presigned URLs:

  K=runs/20261005T210235Z/results/cellranger/cellranger/count/PBMC1K/outs/possorted_genome_bam.bam
  export BAM_URL=$(aws s3 presign s3://<bucket>/$K --expires-in 3600)
  export BAI_URL=$(aws s3 presign s3://<bucket>/$K.bai --expires-in 3600)
  python 05_dhfr_locus.py

(A local BAM path also works for BAM_URL/BAI_URL.)
Output: tables/dhfr_locus_reads.tsv, plus counts of the same genes in both quantifications.
"""

import collections
import os

import numpy as np
import pandas as pd
import pysam

from common import TABLES, load_all

LOCI = {"DHFR": ("chr5", 80626226, 80655002),          # ENSG00000228716, minus strand (GRCh38)
        "lncRNA chr8": ("chr8", 32927913, 33076792)}   # ENSG00000247134, plus strand
GENES = {"ENSG00000228716": "DHFR", "ENSG00000247134": "lncRNA chr8", "ENSG00000113318": "MSH3"}

bam = pysam.AlignmentFile(os.environ["BAM_URL"], index_filename=os.environ["BAI_URL"])
rows = []
for locus, span in LOCI.items():
    mapq, region, nh, counted, low_pos, n = (collections.Counter(), collections.Counter(), collections.Counter(), 0, [], 0)
    for r in bam.fetch(*span):
        if r.is_secondary or r.is_supplementary:
            continue
        n += 1
        tags = dict(r.get_tags())
        mapq[r.mapping_quality] += 1
        region[tags.get("RE", "?")] += 1
        nh[tags.get("NH", 0)] += 1
        counted += bool(tags.get("xf", 0) & 8)  # Cell Ranger: read is the representative of a counted UMI
        if r.mapping_quality < 255:
            low_pos.append(r.reference_start)
    windows = collections.Counter(p // 200 * 200 for p in low_pos).most_common(1)
    rows += [
        (f"{locus}: primary reads in span", n),
        (f"{locus}: unique (MAPQ 255)", mapq[255]),
        (f"{locus}: multi-mapped (MAPQ < 255)", n - mapq[255]),
        (f"{locus}: intronic (RE=N)", region["N"]),
        (f"{locus}: exonic (RE=E)", region["E"]),
        (f"{locus}: reads representing a counted UMI (xf & 8)", counted),
        (f"{locus}: hotspot 200-bp window start (multi-mapped)", windows[0][0] if windows else None),
        (f"{locus}: multi-mapped reads in hotspot window", windows[0][1] if windows else 0),
    ]
    rows += [(f"{locus}: reads with NH={k}", nh[k]) for k in sorted(nh)]
    print(locus, "multi-mapped read positions, 12 bins:", np.histogram(low_pos, bins=12, range=span[1:])[0].tolist())

cr, sa, cr_cells, sa_cells = load_all()
cells = list(cr_cells & sa_cells)
for g, name in GENES.items():
    i = cr.var_names.get_loc(g)
    rows.append((f"{name} UMIs in agreed cells, Cell Ranger", int(cr[cells].X[:, i].sum())))
    for layer in ("spliced", "unspliced", "ambiguous"):
        rows.append((f"{name} UMIs in agreed cells, Simpleaf {layer}", int(sa[cells].layers[layer][:, i].sum())))

out = pd.DataFrame(rows, columns=["metric", "value"])
out.to_csv(TABLES / "dhfr_locus_reads.tsv", sep="\t", index=False)
print(out.to_string(index=False))
