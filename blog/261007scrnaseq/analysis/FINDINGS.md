# 发现：Cell Ranger 与 Simpleaf 为什么保留了不同的细胞

数据：10x PBMC 1k v3，一位健康供者，3′ v3 化学。两条路线都用 nf-core/scrnaseq 4.2.0，处理同一批 FASTQ，基于同一份 GRCh38-2024-A FASTA/GTF。
版本：Cell Ranger 10.0.0；Simpleaf 0.19.5（alevin-fry 0.11.2，`cr-like` 定量）+ QCatch 0.2.12。运行记录见 [../README.md](../README.md) 第 8 节。
每个数字都可以从 `tables/` 里的文件和 `0*_*.py` 脚本重新得到，运行方法见 [README.md](README.md)。

## 一句话结论

两边一致率 99.4%：共同识别 1,221 个，**仅 Simpleaf + QCatch 识别 7 个**，没有仅 Cell Ranger 识别的。这 7 个分歧可以全部追到具体原因：3 个由一个基因位点的定量差异决定，3 个是统计检验的边界细胞，1 个是濒死细胞。它们对细胞类型比例的影响 ≤ 0.17 个百分点。在同一套流程里，只换随机种子或只换计数矩阵，比例变化可达 6–10 个百分点，比细胞识别的影响大一个数量级以上。

## 1. 两个工具做的其实是同一件事

在这份数据上，Cell Ranger 10.0.0 与 QCatch 0.2.12 的细胞识别算法相同，两边的源码都已核对：

1. **OrdMag**：根据 UMI 排在最前面的 barcode 估计细胞数，取 99% 分位数的 1/10 作为初始阈值。
2. **EmptyDrops**：其余 UMI ≥ max(500, 背景最大 UMI + 1) 的 barcode，与环境 RNA 背景做多项分布检验。背景取 UMI 排名第 45,000–90,000 的 barcode，BH 校正后的 FDR 阈值为 0.001。

实现上只有两处不同：Monte Carlo 模拟次数（Cell Ranger 100,000 次，QCatch 10,000 次）和随机数流。Cell Ranger 另有几个额外过滤步骤（高占用 GEM、线粒体阈值、全局最少 UMI），在 3′ v3 基因表达数据上默认都不生效。

**验证**：把 Cell Ranger 的原始矩阵交给 QCatch 的代码（条件 B），得到的 1,221 个细胞与 Cell Ranger 官方结果**逐个 barcode 完全相同**。随机种子 1、2、3、42 和 10,000 / 100,000 次模拟下都是如此。把 Simpleaf 的原始矩阵交给同一份代码（条件 A，100,000 次模拟），结果与 QCatch 官方的 1,228 个也完全相同。

所以分歧**不来自算法**，只来自交给算法的输入：计数和 barcode 范围。

## 2. 输入的两处差异

**barcode 范围**（`tables/raw_universe_summary.tsv`）：

| | barcode 数 | UMI ≥ 100 | UMI ≥ 500 |
| --- | --- | --- | --- |
| Cell Ranger 原始矩阵（≥ 1 个 UMI） | 329,735 | 1,380 | 1,249 |
| Simpleaf 原始矩阵（`--unfiltered-pl`，≥ 10 条 reads） | 72,612 | 1,402 | 1,261 |
| 两边共有 | 71,929 | | |

只在 Cell Ranger 原始矩阵里的 257,806 个 barcode 都是空液滴：中位数 1 个 UMI，最多 82 个，合计只占全部 UMI 的 2.9%。它们的唯一作用是改变 EmptyDrops 背景所取的排名区间。

**定量**（`figures/01_umi_concordance.png`）：在共同识别的 1,221 个细胞上，两边每个细胞的 UMI 数 Spearman 相关 0.9997，Simpleaf 平均多约 5%（log2 比值中位数 0.071）。这里比较的是 Simpleaf 的 `X`，即 spliced + unspliced + ambiguous 之和，已核对完全相等；它与 Cell Ranger 默认包含内含子的计数口径一致。

多出来的那部分并不均匀，集中在少数几个位点（在共同识别细胞里的总计数）：

| 基因 | Cell Ranger | Simpleaf | 说明 |
| --- | --- | --- | --- |
| DHFR（chr5） | 739 | 45,624 | 其中 unspliced 45,532 |
| ENSG00000247134（chr8 lncRNA） | 106 | 23,806 | 其中 unspliced 23,748 |
| MALAT1、线粒体基因 | | 高 5–9% | spliced |
| HLA-C、HLA-DRB5/6、RPSA | | 明显偏低（例如 HLA-DRB6 3,012 → 800） | HLA 有多个高度相似的旁系基因，RPSA 有很多假基因 |

## 3. DHFR 位点（`05_dhfr_locus.py`，`tables/dhfr_locus_reads.tsv`）

直接从 S3 读取 Cell Ranger BAM 中的 DHFR 区域（chr5:80,626,226–80,655,002）：

- 8,988 条 primary reads，其中 **83% 是多重比对**（MAPQ 1 或 3，NH = 2–4），只有 1,519 条唯一比对。
- 多重比对的 reads 几乎全部集中在**一个 200 bp 的内含子窗口**：chr5:80,651,400–80,651,600，有 7,219 条。也就是说，这段序列在基因组中有 2–4 个拷贝。
- Cell Ranger 只统计唯一比对的 reads，最终 DHFR 只有 859 条 reads 对应到计数的 UMI。
- Simpleaf 给 DHFR 记了 45,532 个 unspliced UMI，比 Cell Ranger 在整个区域里看到的 reads 还多约 5 倍。

chr8 那个 lncRNA 的情况不同：Cell Ranger 在它约 150 kb 的区间里一共只有 513 条 reads，几乎全是唯一比对，而 Simpleaf 记了 23,748 个 unspliced UMI。

**解释（推断，未逐条 read 证实）**：这些 reads 在全基因组比对（STAR）里被放到了别的位置，或者因为多重比对被丢弃；Simpleaf 的 splici 索引只含转录本和内含子序列，把它们判给了这两个基因的内含子。Simpleaf 的 RAD 输出不记录每条 read 的基因组坐标，所以没法在 Simpleaf 这一侧直接追踪。

## 4. 七个有争议的细胞

编号与 `figures/01_knee_plots.png` 和 `figures/03_umap_disputed.png` 中的一致。

| # | barcode | 类型（marker） | UMI CR → SA | 线粒体 % | 原因 |
| --- | --- | --- | --- | --- | --- |
| 1 | AGTACTGGTATCGCGC | 血小板（PPBP 30，PF4 3） | 455 → 508 | 5.9 | 越过 500 UMI 门槛：+19 来自 DHFR，+4 来自 lncRNA |
| 5 | GTTCGCTTCGACATTG | 无 marker，受损细胞 | 294 → 577 | 52 | 越过 500 UMI 门槛：+236 来自 DHFR；unspliced 占 53% |
| 6 | TACAACGAGGCCTGAA | 血小板（PPBP 15） | 464 → 532 | 18 | 越过 500 UMI 门槛：+29 来自 DHFR，+8 来自 lncRNA |
| 2 | CTTCAATTCGAGAGCA | cDC（FCER1A、CST3、LYZ） | 1,200 → 1,262 | 4.8 | EmptyDrops 边界：padj 约 5e-3（CR 计数）对比 1e-4（SA 计数） |
| 3 | GGGTCACTCAATCTCT | 单核/树突细胞系（LYZ、LST1、CST3、CD4） | 1,309 → 1,366 | 7.3 | EmptyDrops 边界：padj 约 4e-3 对比 1e-4 |
| 4 | GTGTTAGAGACGGTCA | FCGR3A 单核细胞 | 1,239 → 1,298 | 2.5 | EmptyDrops 边界：padj 约 4e-3 对比 1e-4 |
| 7 | TTGATGGTCGCTAAAC | 无 marker，濒死细胞（67 个基因） | 618 → 716 | 90 | 两边都在边界上：padj 5e-4 到 3e-3 |

这 7 个都是在 EmptyDrops 这一步被加入的，没有一个是 OrdMag 初始识别出来的。也都不像双细胞（Scrublet 得分 ≤ 0.09，在共同识别细胞的分布之内）。

### 交叉实验（`02_cell_calling_decomposition.py`）

用同一份 QCatch 代码，分别替换计数、barcode 范围、模拟次数和随机种子。表格只列 7 个有争议的细胞，"是"表示该条件下被识别为细胞，"未检验"表示 UMI 低于 500：

| 条件 | 计数 | barcode 范围 | 1 | 5 | 6 | 2 | 3 | 4 | 7 | 识别总数 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B | Cell Ranger | Cell Ranger 原始 | 未检验 | 未检验 | 未检验 | 否 | 否 | 否 | 否 | 1,221（= Cell Ranger） |
| G | Cell Ranger 去掉 2 个基因 | Cell Ranger 原始 | 未检验 | 未检验 | 未检验 | 否 | 否 | 否 | 否 | 1,221（= Cell Ranger） |
| C | Cell Ranger | 两边共有 | 未检验 | 未检验 | 未检验 | 否 | **是** | 否 | 否 | 1,221 |
| D | Simpleaf | 两边共有 | 是 | 是 | 是 | 是 | 是 | 是 | 否 | 1,226 |
| E | Simpleaf | Simpleaf 加上 Cell Ranger 的尾部 | 是 | 是 | 是 | 是 | 是 | 否 | 是 | 1,231 |
| F | Simpleaf **去掉 DHFR 和 lncRNA** | Simpleaf 原始 | **未检验** | **未检验** | **未检验** | 是 | 是 | 是 | 否 | 1,225 |
| A | Simpleaf | Simpleaf 原始 | 是 | 是 | 是 | 是 | 是 | 是 | 是 | 1,228（= QCatch，100,000 次模拟） |

从这张表读出：

- **1、5、6 号由 DHFR 和 chr8 lncRNA 这两个位点决定**。去掉这两个基因（F），它们就回到 500 UMI 以下，不再被检验。用 Cell Ranger 计数时，无论 barcode 范围怎么设，它们都没被检验（B、C）。
- **2、3、4 号由其余的定量差异决定，与这两个位点无关**。用 Simpleaf 计数，它们几乎总被识别（A、D、F，padj 1e-4 到 3e-4）；用 Cell Ranger 计数，几乎总不被识别（B、G，padj 4e-3 到 6e-3）。barcode 范围只在边缘起作用：C 里 3 号被识别了，E 里 4 号没被识别。在第 1 至 4 节这些分析的范围内，没有找到"同样的细胞，Simpleaf 计数下离环境 RNA 更远"的具体机制。
- **7 号在所有条件下都在 FDR 阈值附近**，去掉 DHFR（F）或改变 barcode 范围（D）都会让它掉出去。
- **随机性和模拟次数都不改变这 7 个的判定**。A 和 B 在种子 1、2、3、42 下，以及 10,000 / 100,000 次模拟下，对这 7 个的判定完全一致。随机性只让另外 1–2 个边界 barcode 翻转：条件 A 的识别总数在 1,227–1,229 之间。
- 边界本身是模糊的：在各个条件下，至少有一次翻转过的 barcode 还有 7 个（见 `logs/`），并不只是这 7 个。

### 它们像细胞还是像空液滴？（`figures/03_depth_matched.png`）

UMI 数少，表达向量就稀疏，相关系数本身就会偏低。所以每个有争议的细胞都和两组对照比较：下采样到同样 UMI 数的同类型共同识别细胞，以及从环境 RNA 背景里以同样深度模拟的"空液滴"。

- 血小板 1、6 号，以及濒死细胞 7 号：与同类细胞的相关性（0.60–0.72）远高于模拟空液滴（≤ 0.26），明显是细胞。
- 受损细胞 5 号：0.51，介于细胞（下采样后的中位数 0.86）和空液滴（0.22）之间。
- **2、3、4 号（cDC 和单核细胞，约 1,200 UMI）**：与同类细胞的相关性（0.53–0.56）只比模拟空液滴的 95% 分位数（0.49–0.52）高一点。原因是 PBMC 的环境 RNA 本来就以单核细胞的 RNA（LYZ、S100A8/9 等）为主，到这个深度，单核细胞和 DC 很难与环境 RNA 区分开。这正好解释了为什么它们的 p 值卡在阈值附近。

## 5. 对细胞类型比例的影响（`03_disputed_cells.py`，`04_annotation_stability.py`）

同一套聚类加 marker 注释流程。marker 得分低于 0.3 的群标为 "Low quality"，这个群有 36 个细胞，多数线粒体比例超过 90%。最初它在负分里以 0.07 的微小差距被标成了 "Platelet"，后来改正。

| 变化来源 | 比例的最大变化（百分点） |
| --- | --- |
| **细胞识别**（Cell Ranger 细胞集 vs Simpleaf 细胞集；同一计数、同一种子） | **0.16**（标准 QC 后 0.17） |
| 只换随机种子 0–4（同一计数、同一细胞集） | 最大 10.5，中位数 6.8 |
| 只换计数矩阵（同一种子、同一细胞集） | 最大 9.6，中位数 6.5 |

- 标准 QC（≥ 200 个基因、线粒体 ≤ 20%）之后，两种细胞集只差 3 个细胞：2 个 cDC（cDC 比例 1.36% 对 1.53%）和 1 个 FCGR3A 单核细胞。两个血小板因基因数不足 200，两个濒死细胞因线粒体比例过高，都被过滤掉了。
- 换计数矩阵，主要改变的是 NK 和 CD8 T 之间的划分：用 Cell Ranger 计数，NK 为 134 个；用 Simpleaf 计数，NK 为 55 个，5 个种子下都是如此。原因是有一个细胞群的 NK 得分只比 CD8 T 高 0.01。换种子主要改变的是 CD14 和 FCGR3A 单核细胞之间的划分。
- **这些波动幅度来自这套简单的"marker 打分、按群注释"方法**。换用更稳健的注释方法，波动可能会小一些。能下的结论是相对的：在这份数据和这套流程里，细胞识别的差异远小于下游注释的不确定性。

### 一个具体的 marker 结论（`06_marker_conclusion.py`）

用 step 3 的细胞类型标签和 Cell Ranger 计数，对每种细胞集做同样的 Wilcoxon 检验（该类型对比其余所有细胞）：

- **cDC**（16 个对 18 个）：CST3（log2FC 5.2，100% 表达）、FCER1A（6.7）、CLEC10A（5.2）的结论不变。**CD1C 变弱**：cDC 中表达比例从 62.5% 降到 55.6%，校正后 p 从 0.012 升到 0.034。多出的 2 个细胞不表达 CD1C。前 10 名有 3 个基因不同。
- **血小板**：QC 之后两种细胞集里的血小板完全相同（都是 15 个），PPBP 的 log2FC 和校正后 p 也相同，但排名分别是第 13 和第 1。原因是排在前面的基因得分接近并列，参照组变了，排序就跟着变。所以对小细胞群，应比较效应大小，不比较排名。
- **FCGR3A 单核细胞**：FCGR3A、LST1、MS4A7 的 log2FC 基本不变。

## 6. 局限与待做

- 只有一个样本，而且是干净的 PBMC，分歧只有 7 个细胞。这里的结论不能推广到环境 RNA 更多、细胞活性更低或其他组织的样本。下一步可以选 10x PBMC 10k，或者质量较差的样本。
- DHFR 和 chr8 lncRNA 的机制，在 Simpleaf 这一侧是推断。要证实，需要在 Simpleaf 中加入基因组 decoy（如果支持）重新建索引，或者用其他方式拿到 read 级别的比对位置。
- 2、3、4 号为什么在 Simpleaf 计数下更显著，还没有找到单一的原因。
- 用同一种方法重新识别两份原始矩阵（第 1 节）时，差距没有缩小，但这里只用了 Cell Ranger 这一类算法。没有再用 DropletUtils `emptyDrops` 等独立方法做交叉验证。

## 文件索引

| 内容 | 文件 |
| --- | --- |
| 三组 barcode 汇总 | `tables/barcode_sets_summary.tsv` |
| 1,228 个细胞在两份矩阵中的各项指标 | `tables/union_cells_metrics.tsv`，`tables/union_cells_annotated.tsv` |
| 交叉实验 | `tables/cell_calling_conditions_*.tsv`，`tables/cell_calling_barcodes_*.tsv` |
| 有争议细胞的额外 UMI 来自哪些基因 | `tables/disputed_count_differences.tsv` |
| 深度匹配对照 | `tables/disputed_depth_matched.tsv` |
| DHFR / lncRNA 的 read 级证据 | `tables/dhfr_locus_reads.tsv` |
| 细胞类型比例及稳定性 | `tables/celltype_proportions.tsv`，`tables/annotation_stability*.tsv` |
| 图 | `figures/01_knee_plots.png`，`01_umi_concordance.png`，`03_umap_disputed.png`，`03_qc_disputed.png`，`03_depth_matched.png`，`03_marker_dotplot.png` |
