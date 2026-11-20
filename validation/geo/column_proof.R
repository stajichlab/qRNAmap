# column_proof.R: which count columns do the authors' CPM-table columns reproduce, and do the
# replicate assignments fit the data?
#
# Usage (absolute paths):
#   Rscript column_proof.R AUTH_COUNTS_CSV_GZ AUTH_CPM_CSV_GZ ENA_RUNS_TSV COLUMN_MAPPING_TSV \
#       QUANT_COUNTS_ZST OUT_DIR
# Writes OUT_DIR/cpm_column_proof.tsv, OUT_DIR/mapping_totals_vs_ena.tsv,
# OUT_DIR/mapping_spearman_6x6.tsv.
# CPM of a count column = count / column total * 1e6 (all 6,692 genes in the total). Group value =
# log2(mean CPM of the 3 columns of the group). Compared with the CPM-table value on the genes of the
# CPM table where both are finite and not NA. Also the Pearson of the table value with the linear
# group mean CPM (to show the table is on the log2 scale).
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6)
od <- args[6]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
wt <- function(x, f) {
  write.table(x, file.path(od, f), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n==", f, "\n"); print(x, row.names = FALSE)
}
cnt <- read.csv(gzfile(args[1]), check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
cnt <- as.matrix(cnt)
cpmt <- read.csv(gzfile(args[2]), check.names = FALSE, stringsAsFactors = FALSE)
cpm <- t(t(cnt) / colSums(cnt)) * 1e6
groups <- list(
  "2I,4E,9E (AMB only)" = c("2I", "4E", "9E"),
  "3F,5G,8A (control for AMB only)" = c("3F", "5G", "8A"),
  "2F,5C,9D (AMB + LF)" = c("2F", "5C", "9D"),
  "1B,3G,4I (control for AMB + LF)" = c("1B", "3G", "4I"),
  "10B,10D,10F (LF only)" = c("10B", "10D", "10F"),
  "10G,10H,10I (control for LF only)" = c("10G", "10H", "10I"))
rows <- list()
for (tc in setdiff(colnames(cpmt), "systematic_name")) {
  tv <- setNames(cpmt[[tc]], cpmt$systematic_name)
  for (gn in names(groups)) {
    lin <- rowMeans(cpm[names(tv), groups[[gn]], drop = FALSE])
    lg <- log2(lin)
    ok <- !is.na(tv) & is.finite(lg)
    d <- abs(tv[ok] - lg[ok])
    rows[[length(rows) + 1]] <- data.frame(
      cpm_table_column = tc, count_columns = gn, genes = sum(ok),
      median_abs_diff_log2 = median(d), max_abs_diff_log2 = max(d),
      pearson_table_vs_log2_mean_cpm = cor(tv[ok], lg[ok]),
      pearson_table_vs_linear_mean_cpm = cor(tv[ok], lin[ok]))
  }
}
wt(do.call(rbind, rows), "cpm_column_proof.tsv")

# Replicate mapping checks.
mp <- read.delim(args[4], stringsAsFactors = FALSE)
mp <- mp[mp$our_sample != "", ]
ena <- read.delim(args[3], stringsAsFactors = FALSE)
mp$ena_read_count <- ena$read_count[match(mp$srr, ena$run)]
mp$authors_total <- colSums(cnt)[mp$authors_column]
mp$ratio_authors_over_ena <- mp$authors_total / mp$ena_read_count
tot <- mp[, c("authors_column", "our_sample", "srr", "authors_total", "ena_read_count",
              "ratio_authors_over_ena")]
wt(tot, "mapping_totals_vs_ena.tsv")
cat("Spearman totals vs ENA:", cor(tot$authors_total, tot$ena_read_count, method = "spearman"),
    " Pearson:", cor(tot$authors_total, tot$ena_read_count), "\n")
q <- as.matrix(read.delim(pipe(paste("zstd -dc", shQuote(args[5]))), row.names = 1,
                          check.names = FALSE, quote = ""))
g <- intersect(rownames(q), rownames(cnt))
cat("common genes:", length(g), "\n")
sm <- sapply(colnames(q), function(s)
  sapply(mp$authors_column, function(a) cor(cnt[g, a], q[g, s], method = "spearman")))
rownames(sm) <- paste0(mp$authors_column, "=", mp$our_sample)
out <- data.frame(authors_column = rownames(sm), sm, check.names = FALSE)
wt(out, "mapping_spearman_6x6.tsv")
best <- colnames(sm)[apply(sm, 1, which.max)]
cat("best match per authors' column:", paste(rownames(sm), best, sep = " -> "), "\n")
cond <- sub("_.*", "", mp$our_sample)
dg <- diag(sm[, mp$our_sample])
same <- outer(cond, cond, "==") & !diag(6)
smo <- sm[, mp$our_sample]
cat("diagonal:", range(dg), "\nsame condition off-diagonal:", range(smo[same]),
    "\nacross conditions:", range(smo[outer(cond, cond, "!=")]), "\n")
cat("per row: diagonal minus best same-condition off-diagonal:",
    sapply(1:6, function(i) dg[i] - max(smo[i, same[i, ]])), "\n")
