# compare_modes.R QUANT_DIR ALIGN_DIR OUT_DIR  (copy of validation/compare_modes.R, same logic)
# Compares the gene count matrices and the DESeq2 results of the two qrnamap modes.
# Inputs: DIR/counts.tsv.zst (qrnamap matrix) and DIR/deseq2.results.csv (qrnamap run --contrast).
# Spearman: R stats::cor(method = "spearman") (ties get average ranks).
# Writes OUT_DIR/mode_agreement.tsv and OUT_DIR/de_overlap.tsv (small tables) and prints them.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
qd <- args[1]; ad <- args[2]; od <- args[3]
dir.create(od, showWarnings = FALSE, recursive = TRUE)

read_counts <- function(d) {
  f <- file.path(d, "counts.tsv.zst")
  con <- pipe(paste("zstd -dc", shQuote(f)))
  m <- read.delim(con, row.names = 1, check.names = FALSE, quote = "")
  as.matrix(m)
}
q <- read_counts(qd)
a <- read_counts(ad)
cat("rows quant", nrow(q), "align", nrow(a), "\n")
cat("identical row names (same order):", identical(rownames(q), rownames(a)), "\n")
cat("identical column names:", identical(colnames(q), colnames(a)), "\n")
stopifnot(identical(rownames(q), rownames(a)), identical(colnames(q), colnames(a)))

rows <- list()
for (s in colnames(q)) {
  x <- q[, s]; y <- a[, s]
  both0 <- x == 0 & y == 0
  rows[[s]] <- data.frame(
    sample = s, genes = length(x),
    spearman_all_genes = cor(x, y, method = "spearman"),
    genes_nonzero_either = sum(!both0),
    spearman_nonzero_either = cor(x[!both0], y[!both0], method = "spearman"),
    pearson_log2p1 = cor(log2(x + 1), log2(y + 1)),
    total_quant = sum(x), total_align = sum(y),
    zero_quant = sum(x == 0), zero_align = sum(y == 0))
}
agree <- do.call(rbind, rows)
pooled <- data.frame(
  sample = "ALL_pooled", genes = length(q),
  spearman_all_genes = cor(as.vector(q), as.vector(a), method = "spearman"),
  genes_nonzero_either = sum(!(q == 0 & a == 0)),
  spearman_nonzero_either = cor(as.vector(q)[!(q == 0 & a == 0)], as.vector(a)[!(q == 0 & a == 0)], method = "spearman"),
  pearson_log2p1 = cor(log2(as.vector(q) + 1), log2(as.vector(a) + 1)),
  total_quant = sum(q), total_align = sum(a), zero_quant = sum(q == 0), zero_align = sum(a == 0))
agree <- rbind(agree, pooled)
# Per gene across the 6 samples (genes with a nonzero count in both modes and some variation).
pg <- sapply(seq_len(nrow(q)), function(i) {
  x <- q[i, ]; y <- a[i, ]
  if (sd(x) == 0 || sd(y) == 0) NA_real_ else cor(x, y, method = "spearman")
})
cat("per-gene Spearman across samples: genes with a value", sum(!is.na(pg)),
    "median", median(pg, na.rm = TRUE), "\n")
write.table(agree, file.path(od, "mode_agreement.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(agree, digits = 6)

read_de <- function(d) read.csv(file.path(d, "deseq2.results.csv"), check.names = FALSE)
dq <- read_de(qd); da <- read_de(ad)
sig <- function(r) r$gene[!is.na(r$padj) & r$padj < 0.05]
sq <- sig(dq); sa <- sig(da)
inter <- intersect(sq, sa); uni <- union(sq, sa)
m <- merge(dq[, c("gene", "log2FoldChange", "padj")], da[, c("gene", "log2FoldChange", "padj")],
           by = "gene", suffixes = c("_quant", "_align"))
ii <- m$gene %in% inter
same_sign <- sum(sign(m$log2FoldChange_quant[ii]) == sign(m$log2FoldChange_align[ii]))
lfc_ok <- !is.na(m$log2FoldChange_quant) & !is.na(m$log2FoldChange_align)
ov <- data.frame(
  metric = c("genes_tested_quant", "genes_tested_align", "padj_na_quant", "padj_na_align",
             "sig_quant", "sig_align", "sig_both", "sig_only_quant", "sig_only_align",
             "jaccard", "sig_both_same_lfc_sign", "spearman_lfc_all_genes_with_lfc"),
  value = c(nrow(dq), nrow(da), sum(is.na(dq$padj)), sum(is.na(da$padj)),
            length(sq), length(sa), length(inter), length(setdiff(sq, sa)), length(setdiff(sa, sq)),
            length(inter) / length(uni), same_sign,
            cor(m$log2FoldChange_quant[lfc_ok], m$log2FoldChange_align[lfc_ok], method = "spearman")))
write.table(ov, file.path(od, "de_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(ov, digits = 6)
