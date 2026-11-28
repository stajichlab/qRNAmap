#!/usr/bin/env Rscript
# Check the precomputed prior histograms (deseq2_prior_table.R) against DESeq2: for every dump directory
# (deseq2_dump.R) whose dispPriorVar came from the m - p <= 3 simulation, recompute dispPriorVar from genes.tsv and the
# table (KL divergence of the residual histogram, loess(span = 0.2), argmin on 1000 points, floor 0.25), as
# estimateDispersionsPriorVar does, and require the value DESeq2 reported.
#   Rscript check_prior_table.R TABLE.tsv DUMP_DIR...
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("usage: check_prior_table.R TABLE.tsv DUMP_DIR...")
tab <- read.delim(a[1])
brks <- -20:20 / 2
grid <- seq(from = 0, to = 8, length = 200)
stopifnot(isTRUE(all.equal(tab$x[tab$df == 1], grid, tolerance = 0)))
bad <- 0; checked <- 0
for (d in a[-1]) {
  run <- read.delim(file.path(d, "run.tsv"), header = FALSE, col.names = c("key", "value"), quote = "")
  val <- setNames(run$value, run$key)
  if (val[["dispPriorVar_path"]] != "simulation") next
  df <- as.integer(val[["m"]]) - as.integer(val[["p"]])
  g <- read.delim(file.path(d, "genes.tsv"), check.names = FALSE)
  g <- g[!g$allZero, ]
  above <- g$dispGeneEst >= 1e-8 * 100
  obs <- (log(g$dispGeneEst) - log(g$dispFit))[above]
  obs <- obs[obs > min(brks) & obs < max(brks)]
  od <- hist(obs, breaks = brks, plot = FALSE)$density
  t <- tab[tab$df == df, ]
  kl <- vapply(seq_len(nrow(t)), function(i) {
    rd <- as.numeric(t[i, paste0("c", 1:40)]) / (t$n[i] * 0.5)
    z <- c(od, rd); small <- min(z[z > 0])
    sum(od * (log(od + small) - log(rd + small)))
  }, 0)
  lofit <- loess(kl ~ grid, span = 0.2)
  fine <- seq(from = 0, to = 8, length = 1000)
  got <- pmax(fine[which.min(predict(lofit, fine))], 0.25)
  want <- as.numeric(val[["dispPriorVar"]])
  checked <- checked + 1
  ok <- abs(got - want) <= 1e-12 * max(1, abs(want))
  if (!ok) bad <- bad + 1
  cat(sprintf("%s %s: m - p = %d, DESeq2 %.17g, from the table %.17g\n", if (ok) "ok" else "FAIL", d, df, want, got))
}
if (checked == 0) stop("no dump directory used the m - p <= 3 simulation")
if (bad > 0) quit(status = 1)
