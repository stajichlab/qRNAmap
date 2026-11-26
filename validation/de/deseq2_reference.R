#!/usr/bin/env Rscript
# DESeq2 reference results for a count matrix, for comparing the page's DE (site/analyze/de-stats.js) and,
# later, a Rust port against R DESeq2 (pinned in pixi.toml).
#   Rscript deseq2_reference.R COUNTS.tsv COLDATA.tsv DESIGN FACTOR LEVEL REFERENCE OUT_PREFIX
# COUNTS: gene column then one integer column per sample. COLDATA: column `sample` plus covariates.
# DESIGN: an R formula such as "~ condition" or "~ strain + condition".
# Writes OUT_PREFIX.results.tsv (gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj; DESeq2 defaults:
# Wald test, independent filtering, Cook's cutoff) and OUT_PREFIX.sizefactors.tsv.
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 7) stop("usage: deseq2_reference.R COUNTS COLDATA DESIGN FACTOR LEVEL REFERENCE OUT_PREFIX")
counts <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
storage.mode(counts) <- "integer" # kallisto est_counts summed per gene: rounded, as tximport does for DESeq2
cd <- read.delim(a[2], check.names = FALSE, stringsAsFactors = TRUE)
rownames(cd) <- cd$sample
cd <- cd[colnames(counts), , drop = FALSE]
cd[[a[4]]] <- relevel(factor(cd[[a[4]]]), ref = a[6])
dds <- DESeqDataSetFromMatrix(counts, cd, as.formula(a[3]))
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c(a[4], a[5], a[6]))
out <- data.frame(gene = rownames(res), as.data.frame(res)[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")])
write.table(out, paste0(a[7], ".results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(data.frame(sample = colnames(dds), sizeFactor = sizeFactors(dds)), paste0(a[7], ".sizefactors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("DESeq2 %s, %s: %d genes, padj < 0.05: %d\n", packageVersion("DESeq2"), a[3], nrow(out), sum(out$padj < 0.05, na.rm = TRUE)))
