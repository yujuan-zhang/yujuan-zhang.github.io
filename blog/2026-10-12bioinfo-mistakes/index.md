# Small Mistakes, Big Consequences: Notes on Common Errors in Bioinformatics

Most bioinformatics errors are not exotic. They are an index that starts at the wrong number, a genome build nobody wrote down, or a stray space in a shell command. They are easy to make, hard to spot, and often fail silently: the pipeline finishes and the output looks reasonable.

## 1. Counting Is Harder Than It Looks

**Coordinate systems differ between file formats.** BED intervals are 0-based and half-open, while GTF and GFF are 1-based and closed. The first base of a chromosome is `0` in BED and `1` in GTF. Converting between them without adjusting the start position shifts every feature by one base, which is enough to change a codon or miss a splice site.

**Programming languages use different index origins.** Python counts from 0; R counts from 1. Code ported from one language to the other can run without errors and still return the neighbouring element.

**Chromosome names should not be sorted as plain text.** A plain text sort places `chr10` before `chr2`, breaking the usual numeric order. Use a natural sort to obtain `chr1`, `chr2`, …, `chr10`; for interval files such as BED, also sort by start coordinate within each chromosome.

> **Good practice:** When a coordinate crosses a format or language boundary, check one known feature by hand, such as the start of a well-studied gene, before trusting the conversion.

## 2. The Data Do Not Match What You Think

**Genome builds must stay consistent.** Reads aligned to hg19 but viewed against hg38 annotation appear to fall outside exons. Every reference file, index and annotation must come from the same assembly and, ideally, the same release.

**Strand matters.** A feature on the minus strand must be reverse-complemented before its sequence is used. Forgetting this produces a valid-looking but wrong sequence.

**Define how to handle reads that match multiple genes before counting.** Genes can overlap or be nested. Without a rule for shared regions, features may be missed or reads counted twice.

**Line endings are a hidden trap.** Files saved on Windows end lines with `\r\n`; Unix tools expect `\n`. A missing newline at the end of the last line can also cause a simple loop to skip it. `dos2unix` and a quick `tail` check help.

> **Good practice:** Record the genome build and annotation version in file names and in a project README.

## 3. Shell Errors

**Deleting or overwriting by accident.**

- A stray space turns `rm *.fasta` into `rm * .fasta`, which deletes everything in the folder.
- `>` replaces a file's contents, while `>>` appends. Redirecting with `>` into `~/.bashrc` erases the existing configuration.
- `rm *` in the wrong directory is not reversible.

**Letting the shell interpret characters meant for a program.**

- `awk {print $2}` fails because the shell expands `$2` and the braces; write `awk '{print $2}'`.
- `grep > file.fasta` redirects output instead of searching for headers; write `grep '>' file.fasta`.
- `find . -name *.txt` may expand the wildcard before `find` sees it; quote the pattern.

**Forgetting a default.**

- `cut` splits on tabs unless told otherwise; a CSV needs `-d,`.
- `mkdir` needs `-p` to create nested directories.
- `rm` needs `-r` for directories.
- `cat` on a multi-gigabyte FASTQ floods the terminal; `head` or `less -S` is safer.
- `chmod 777` grants everyone full access; `chmod u+x` is usually all that is needed.

> **Good practice:** Before deleting with a wildcard, run `ls` with the same pattern first. In Bash, `set -o noclobber` stops `>` from overwriting an existing file.

## 4. Process Mistakes Outlast Code Mistakes

**Use established tools rather than reinventing the wheel.** To read FASTQ, GenBank or BLAST output, prefer tested libraries such as Biopython, pysam or Bioconductor packages: they save effort and handle edge cases that custom parsers often miss. For single-cell RNA-seq, use Scanpy or Seurat for quality control, dimensionality reduction and clustering, and spend the time saved on the research question.

**Test regular expressions** against examples you know should and should not match.

**Track versions of data and code.** Analyses that silently use an outdated input are hard to reproduce. Git, Git LFS or a workflow manager makes the inputs explicit.

**Document decisions and use appropriate statistics.** A result that cannot be explained, or that depends on an unsuitable test or poor randomization, is not yet a finding.

## 5. Common Mistakes in Single-Cell Analysis

**Treating cells as independent replicates in differential expression.** When comparing conditions such as healthy and disease, each with several samples, a common approach is to pool all cells per condition and test at the cell level. Cells from the same sample are not independent, and thousands of cells make p-values artificially small, producing many false positives. **Pseudobulk analysis is more reliable:** aggregate counts per sample, then compare samples.

**Over-interpreting UMAP.** UMAP is a non-linear dimensionality reduction method; distances between points or clusters have no reliable quantitative meaning. Use it to visualize structure, not to judge how similar cell populations are.

**Assuming batch correction is safe.** Batch correction and data integration rely on assumptions, and over-correction can remove genuine biological signal. A benchmark of seven widely used methods found that many introduce measurable artifacts; only Harmony performed consistently well, and it was the only method the authors recommended for scRNA-seq batch correction.[Antonsson & Melsted, 2024](https://www.biorxiv.org/content/10.1101/2024.03.19.585562v1.full)

> **Good practice:** Test differential expression at the sample level, treat UMAP as visualization only, and confirm that known biological differences survive batch correction.

## A Short Checklist

| Before you... | Check |
| --- | --- |
| Convert coordinates | 0-based or 1-based? Half-open or closed? |
| Combine files | Same genome build and annotation release? Same chromosome naming (`chr1` vs `1`)? |
| Extract sequences | Which strand? |
| Run `rm`, `>` or `mv` | Correct directory? Pattern confirmed with `ls`? |
| Pass code to `awk`, `grep` or `find` | Quoted? |
| Compare conditions | Analysed at the sample level (pseudobulk)? |
| Share results | Inputs, versions and parameters documented? |

None of these mistakes requires advanced knowledge to avoid. They require a moment of checking, which is cheaper than rerunning an analysis.

---

*The main content of this post is based on my reading of Tommy Tang's article "[The Most Common Stupid Mistakes In Bioinformatics](https://divingintogeneticsandgenomics.com/post/stupid-mistakes-for-bioinformatics/)" on his blog, [Chatomics](https://divingintogeneticsandgenomics.com/). Please read the original for his full list and additional resources.*
