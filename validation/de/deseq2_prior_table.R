#!/usr/bin/env Rscript
# The simulated histograms of DESeq2's prior variance for m - p <= 3 (estimateDispersionsPriorVar, DESeq2 1.50.2).
#   Rscript deseq2_prior_table.R OUT.tsv
# For m - p <= 3 and m > p, DESeq2 sets set.seed(2) and, for x on seq(0, 8, length = 200) in order, draws
# log(rchisq(1e4, m - p)) + rnorm(1e4, 0, sqrt(x)) - log(m - p), keeps the values strictly inside (-10, 10) and bins them
# with hist(breaks = -20:20/2). Nothing before the draws uses the random stream, so the histograms depend only on m - p:
# this script repeats that loop (same code, same order) for m - p = 1, 2, 3.
# OUT.tsv: df, i (1..200), x, n (values kept), c1..c40 (counts per bin, left to right). Density = c / (n * 0.5), which is
# what hist() computes. check_prior_table.R checks the table against DESeq2's dispPriorVar on the fixtures.
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 1) stop("usage: deseq2_prior_table.R OUT.tsv")
brks <- -20:20 / 2
grid <- seq(from = 0, to = 8, length = 200)
rows <- list()
for (df in 1:3) {
  set.seed(2)
  for (i in seq_along(grid)) {
    x <- grid[i]
    randDist <- log(rchisq(10000, df = df)) + rnorm(10000, 0, sqrt(x)) - log(df)
    randDist <- randDist[randDist > min(brks) & randDist < max(brks)]
    h <- hist(randDist, breaks = brks, plot = FALSE)
    rows[[length(rows) + 1]] <- c(df, i, sprintf("%.17g", x), length(randDist), h$counts)
  }
}
tab <- do.call(rbind, rows)
colnames(tab) <- c("df", "i", "x", "n", paste0("c", 1:40))
write.table(tab, a[1], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d rows (R %s)\n", a[1], nrow(tab), getRversion()))
