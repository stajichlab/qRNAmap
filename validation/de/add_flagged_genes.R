#!/usr/bin/env Rscript
# Add to a fixture's counts the genes that its full-data DESeq2 run sent to the GLM optim fallback (optimRow) or flagged
# by Cook's distance (p-value NA with a statistic), at most 10 of each, keeping the full matrix's gene order. A random
# subset rarely holds these genes, and the port's stage tests need them.
#   Rscript add_flagged_genes.R FULL_DUMP_DIR FULL_COUNTS.tsv FIXTURE_COUNTS.tsv   (rewrites FIXTURE_COUNTS.tsv)
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) stop("usage: add_flagged_genes.R FULL_DUMP_DIR FULL_COUNTS.tsv FIXTURE_COUNTS.tsv")
g <- read.delim(file.path(a[1], "genes.tsv"), check.names = FALSE)
full <- read.delim(a[2], check.names = FALSE)
fix <- read.delim(a[3], check.names = FALSE)
optim <- head(g$gene[g$optimRow], 10)
cooks <- head(g$gene[!g$allZero & is.na(g$res_pvalue) & !is.na(g$res_stat)], 10)
add <- setdiff(c(optim, cooks), fix[[1]])
out <- full[full[[1]] %in% c(fix[[1]], add), names(fix)]
write.table(out, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: +%d genes (optim %d, Cook's %d in the full run)\n", a[3], length(add), length(optim), length(cooks)))
