# compare_geo.R: compare qRNAmap AMB results with the authors' GSE80357 tables.
#
# Usage (all absolute paths):
#   Rscript compare_geo.R AUTH_COUNTS QUANT_COUNTS ALIGN_COUNTS AUTH_DE QUANT_DE ALIGN_DE \
#       AUTHDESEQ_DE QUANT_RERUN_DE AUTH_CPM OUT_DIR
#   AUTH_COUNTS     authors_counts_amb.tsv (prepare_geo.py): gene + our 6 sample names
#   QUANT_COUNTS, ALIGN_COUNTS   qrnamap counts.tsv.zst of run_quant and run_align
#   AUTH_DE         GSE80357_log_fold_change_and_FDR_values.csv.gz (authors' DE table)
#   QUANT_DE, ALIGN_DE           deseq2.results.csv of run_quant and run_align
#   AUTHDESEQ_DE    results.csv of `qrnamap de` on AUTH_COUNTS
#   QUANT_RERUN_DE  results.csv of `qrnamap de` on QUANT_COUNTS with the binary of this task
#   AUTH_CPM        GSE80357_average_cpm_per_treatment_group.csv.gz (only for a column check)
#   OUT_DIR         small result tables (TSV)
# Definitions:
#   Spearman: cor(method = "spearman") (ties get average ranks). Significant: FDR or padj < 0.05;
#   a gene with NA FDR/padj is not significant. Common genes: present in both tables.
#   log2FC comparisons use common genes with a non-NA log2FC in both tables.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 10)
names(args) <- c("auth_counts", "quant_counts", "align_counts", "auth_de", "quant_de",
                 "align_de", "authdeseq_de", "quant_rerun_de", "auth_cpm", "out")
od <- args[["out"]]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
wt <- function(x, f) {
  write.table(x, file.path(od, f), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n==", f, "\n"); print(x, row.names = FALSE)
}

read_zst <- function(f) {
  m <- read.delim(pipe(paste("zstd -dc", shQuote(f))), row.names = 1, check.names = FALSE,
                  quote = "", colClasses = c("character", rep("numeric", 6)))
  as.matrix(m)
}
au <- as.matrix(read.delim(args[["auth_counts"]], row.names = 1, check.names = FALSE,
                           quote = "", colClasses = c("character", rep("numeric", 6))))
q <- read_zst(args[["quant_counts"]])
a <- read_zst(args[["align_counts"]])
stopifnot(identical(colnames(au), colnames(q)), identical(colnames(q), colnames(a)),
          identical(rownames(q), rownames(a)))
g <- intersect(rownames(q), rownames(au))
cat("common genes:", length(g), "\n")

# ---- 2. counts agreement on the common genes ----
cnt_row <- function(lbl, s, x, y, ref) {
  # x = reference side ("theirs" or align), y = compared side ("ours" or quant)
  keep <- ref >= 10
  r <- (y[keep] + 1) / (x[keep] + 1)
  qs <- quantile(r, c(0.25, 0.5, 0.75), names = FALSE)
  data.frame(comparison = lbl, sample = s, genes = length(x),
             spearman_raw = cor(x, y, method = "spearman"),
             spearman_log2p1 = cor(log2(x + 1), log2(y + 1), method = "spearman"),
             pearson_log2p1 = cor(log2(x + 1), log2(y + 1)),
             total_x = sum(x), total_y = sum(y), total_ratio_y_over_x = sum(y) / sum(x),
             genes_authors_ge10 = sum(keep), ratio_q25 = qs[1], ratio_median = qs[2],
             ratio_q75 = qs[3], ratio_iqr = qs[3] - qs[1])
}
rows <- list()
for (s in colnames(q)) {
  A <- au[g, s]; Q <- q[g, s]; L <- a[g, s]
  rows[[length(rows) + 1]] <- cnt_row("authors(x)_vs_quant(y)", s, A, Q, A)
  rows[[length(rows) + 1]] <- cnt_row("authors(x)_vs_align(y)", s, A, L, A)
  rows[[length(rows) + 1]] <- cnt_row("align(x)_vs_quant(y)", s, L, Q, A)
}
A <- as.vector(au[g, ]); Q <- as.vector(q[g, ]); L <- as.vector(a[g, ])
rows[[length(rows) + 1]] <- cnt_row("authors(x)_vs_quant(y)", "pooled", A, Q, A)
rows[[length(rows) + 1]] <- cnt_row("authors(x)_vs_align(y)", "pooled", A, L, A)
rows[[length(rows) + 1]] <- cnt_row("align(x)_vs_quant(y)", "pooled", L, Q, A)
cnt <- do.call(rbind, rows)
cnt <- cnt[order(cnt$comparison), ]
wt(cnt, "counts_agreement.tsv")
# Reproduce the earlier quant-vs-align Spearman on all 6,477 genes (check).
chk <- sapply(colnames(q), function(s) cor(q[, s], a[, s], method = "spearman"))
cat("\nquant vs align Spearman, all", nrow(q), "genes (check vs validation/amb):\n")
print(round(chk, 6))

# ---- 3./4. DE tables ----
ade <- read.csv(gzfile(args[["auth_de"]]), check.names = FALSE, stringsAsFactors = FALSE)
cat("\nauthors' DE header:", paste(colnames(ade), collapse = ", "), "\n")
cat("authors' DE rows:", nrow(ade), " duplicated genes:", sum(duplicated(ade$genes)), "\n")
for (cn in grep("^AMB_only_", colnames(ade), value = TRUE))
  cat(cn, ": class", class(ade[[cn]]), ", NA", sum(is.na(ade[[cn]])), "\n")
auth <- data.frame(gene = ade$genes, lfc = ade$AMB_only_logFC, fdr = ade$AMB_only_FDR,
                   p = ade$AMB_only_PValue)
read_dq <- function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  data.frame(gene = d$gene, lfc = d$log2FoldChange, fdr = d$padj, p = d$pvalue)
}
dq <- read_dq(args[["quant_de"]]); dl <- read_dq(args[["align_de"]])
dd <- read_dq(args[["authdeseq_de"]]); dr <- read_dq(args[["quant_rerun_de"]])

# Check: the rerun of qrnamap de on the quant counts with this task's binary.
m <- merge(dq, dr, by = "gene")
cat("\nrerun check: genes", nrow(dq), nrow(dr), "merged", nrow(m),
    " max |lfc diff|", max(abs(m$lfc.x - m$lfc.y), na.rm = TRUE),
    " padj NA pattern identical", identical(is.na(m$fdr.x), is.na(m$fdr.y)),
    " max |padj diff|", max(abs(m$fdr.x - m$fdr.y), na.rm = TRUE), "\n")

de_cmp <- function(lbl, x, y) {
  m <- merge(x, y, by = "gene", suffixes = c(".x", ".y"))
  ok <- !is.na(m$lfc.x) & !is.na(m$lfc.y)
  sx <- !is.na(m$fdr.x) & m$fdr.x < 0.05
  sy <- !is.na(m$fdr.y) & m$fdr.y < 0.05
  both <- sx & sy
  either <- (sx | sy) & ok
  data.frame(comparison = lbl, genes_x = nrow(x), genes_y = nrow(y), common = nrow(m),
             common_lfc_both = sum(ok),
             spearman_lfc = cor(m$lfc.x[ok], m$lfc.y[ok], method = "spearman"),
             pearson_lfc = cor(m$lfc.x[ok], m$lfc.y[ok]),
             na_fdr_x = sum(is.na(x$fdr)), na_fdr_y = sum(is.na(y$fdr)),
             sig_x_all = sum(!is.na(x$fdr) & x$fdr < 0.05),
             sig_y_all = sum(!is.na(y$fdr) & y$fdr < 0.05),
             sig_x_common = sum(sx), sig_y_common = sum(sy), sig_both = sum(both),
             x_only = sum(sx & !sy), y_only = sum(sy & !sx),
             jaccard = sum(both) / sum(sx | sy),
             both_same_sign = sum(both & ok & sign(m$lfc.x) == sign(m$lfc.y)),
             either_sig_lfc_both = sum(either),
             either_same_sign = sum(either & sign(m$lfc.x) == sign(m$lfc.y)))
}
de <- rbind(
  de_cmp("authors_published(x)_vs_deseq2_quant(y)", auth, dq),
  de_cmp("authors_published(x)_vs_deseq2_align(y)", auth, dl),
  de_cmp("deseq2_align(x)_vs_deseq2_quant(y)", dl, dq),
  de_cmp("i:authors_published(x)_vs_deseq2_authors_counts(y)", auth, dd),
  de_cmp("ii:deseq2_authors_counts(x)_vs_deseq2_quant(y)", dd, dq),
  de_cmp("ii:deseq2_authors_counts(x)_vs_deseq2_align(y)", dd, dl),
  # The same comparisons restricted to the genes of the authors' DE table (their filter), so
  # that every row below uses about the same gene set as the rows with the published table.
  de_cmp("on_authors_DE_genes:deseq2_align(x)_vs_deseq2_quant(y)",
         dl[dl$gene %in% auth$gene, ], dq[dq$gene %in% auth$gene, ]),
  de_cmp("on_authors_DE_genes:ii:deseq2_authors_counts(x)_vs_deseq2_quant(y)",
         dd[dd$gene %in% auth$gene, ], dq[dq$gene %in% auth$gene, ]),
  de_cmp("on_authors_DE_genes:ii:deseq2_authors_counts(x)_vs_deseq2_align(y)",
         dd[dd$gene %in% auth$gene, ], dl[dl$gene %in% auth$gene, ]))
wt(de, "de_agreement.tsv")

# ---- sign convention ----
cpm <- t(t(au) / colSums(au)) * 1e6
amb <- grep("^AmpB_", colnames(au)); ctl <- grep("^Control_", colnames(au))
lr <- log2((rowMeans(cpm[, amb]) + 0.5) / (rowMeans(cpm[, ctl]) + 0.5))
sc <- auth[auth$gene %in% names(lr), ]
sc$log2_cpm_ratio_ampb_over_control <- lr[sc$gene]
s1 <- sc$fdr < 0.05
sgn <- data.frame(
  check = c("authors' counts: log2((mean CPM AmpB+0.5)/(mean CPM Control+0.5)) vs AMB_only_logFC",
            "our DESeq2 quant (AmpB vs Control) log2FC vs AMB_only_logFC"),
  genes = c(nrow(sc), NA), spearman = c(cor(sc$lfc, sc$log2_cpm_ratio_ampb_over_control,
                                             method = "spearman"), NA),
  authors_fdr_lt_0.05 = c(sum(s1), NA),
  same_sign_among_them = c(sum(sign(sc$lfc[s1]) == sign(sc$log2_cpm_ratio_ampb_over_control[s1])),
                           NA))
mq <- merge(auth, dq, by = "gene", suffixes = c(".auth", ".quant"))
okq <- !is.na(mq$lfc.quant) & mq$fdr.auth < 0.05
sgn$genes[2] <- sum(!is.na(mq$lfc.quant))
sgn$spearman[2] <- cor(mq$lfc.auth, mq$lfc.quant, method = "spearman", use = "complete.obs")
sgn$authors_fdr_lt_0.05[2] <- sum(okq)
sgn$same_sign_among_them[2] <- sum(sign(mq$lfc.auth[okq]) == sign(mq$lfc.quant[okq]))
wt(sgn, "sign_convention.tsv")
top <- head(mq[order(mq$fdr.auth), c("gene", "lfc.auth", "fdr.auth", "lfc.quant", "fdr.quant")], 10)
top$log2_cpm_ratio_authors_counts <- lr[top$gene]
wt(top, "sign_convention_top10.tsv")

# ---- observations on the authors' tables (not a comparison with qRNAmap) ----
cp <- read.csv(gzfile(args[["auth_cpm"]]), check.names = FALSE, stringsAsFactors = FALSE)
cat("\nauthors' CPM header:", paste(colnames(cp), collapse = ", "), " rows", nrow(cp), "\n")
obs <- data.frame(
  check = c("DE table: rows with NA in LF_only_logFC",
            "CPM table: rows with NA in LF_only",
            "DE table: Pearson AMB_only_logFC vs LF_only_logFC (complete rows)",
            "DE table: Pearson AMB_only_logFC vs AMB_and_LF_logFC",
            "DE table: median |AMB_only_logFC - LF_only_logFC| (complete rows)",
            "DE table: max |AMB_only_logFC - LF_only_logFC| (complete rows)",
            "CPM table: max |AMB_only - LF_only| (complete rows)",
            "CPM table: max |control_for_AMB_only - control_for_LF_only| (complete rows)",
            "CPM table: max |AMB_only - AMB_LF|",
            "DE table genes also in counts table",
            "counts table genes in DE table"),
  value = c(sum(is.na(ade$LF_only_logFC)), sum(is.na(cp$LF_only)),
            cor(ade$AMB_only_logFC, ade$LF_only_logFC, use = "complete.obs"),
            cor(ade$AMB_only_logFC, ade$AMB_and_LF_logFC),
            median(abs(ade$AMB_only_logFC - ade$LF_only_logFC), na.rm = TRUE),
            max(abs(ade$AMB_only_logFC - ade$LF_only_logFC), na.rm = TRUE),
            max(abs(cp$AMB_only - cp$LF_only), na.rm = TRUE),
            max(abs(cp$control_for_AMB_only - cp$control_for_LF_only), na.rm = TRUE),
            max(abs(cp$AMB_only - cp$AMB_LF)),
            sum(ade$genes %in% rownames(au)), sum(rownames(au) %in% ade$genes)))
wt(obs, "authors_table_checks.tsv")
# Authors' counts of the genes absent from their DE table.
nd <- setdiff(rownames(au), ade$genes)
cat("\ngenes in authors' counts but not in their DE table:", length(nd),
    "; all-zero in the 6 AMB columns:", sum(rowSums(au[nd, ]) == 0),
    "; max row sum over the 6 AMB columns:", max(rowSums(au[nd, ])), "\n")
cat("genes in authors' DE table: min row sum over the 6 AMB columns:",
    min(rowSums(au[intersect(ade$genes, rownames(au)), ])), "\n")
writeLines(capture.output(sessionInfo()), file.path(od, "compare_sessioninfo.txt"))
