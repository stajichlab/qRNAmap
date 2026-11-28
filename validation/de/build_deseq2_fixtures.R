#!/usr/bin/env Rscript
# Inputs of the DESeq2 port's fixtures (and of the full-data reference runs) from three public RNA-seq datasets.
#   Rscript build_deseq2_fixtures.R RAW_DIR OUT_DIR
# RAW_DIR (make_deseq2_fixtures.sh fetches it):
#   yeast/*.gbgout, exclude.lst   Schurch et al. 2016 (RNA 22:839), S. cerevisiae WT vs snf2, 48 replicates each, htseq
#                                 counts and the authors' bad-replicate list (github.com/bartongroup/profDGE48, MIT)
#   pasilla/pasilla_gene_counts.tsv, pasilla/pasilla_sample_annotation.csv
#                                 Brooks et al. 2011, D. melanogaster pasilla knock-down (Bioconductor pasilla, LGPL)
#   afum/counts.full.tsv, afum/coldata.tsv
#                                 A. fumigatus AF293/CEA10 Azole vs GMM, kallisto counts streamed by the qRNAmap page
#                                 (validation/de/README.md, 2026-11-26)
# OUT_DIR/fixtures/<name>/ and OUT_DIR/full/<name>/ get counts.tsv, coldata.tsv and spec.tsv (design, contrast, ref,
# what the case is for). Fixtures are a seeded subset of genes (300 drawn at random, plus all-zero genes, genes zero in
# one contrasted group, and genes zero in both contrasted levels when present), so DESeq2 runs on exactly the fixture.
# Sample choices and synthetic covariates are seeded; everything here is deterministic.
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 2) stop("usage: build_deseq2_fixtures.R RAW_DIR OUT_DIR")
raw <- a[1]; outdir <- a[2]

## Yeast: one column per replicate; drop htseq's summary rows and the authors' bad replicates.
files <- sort(list.files(file.path(raw, "yeast"), pattern = "\\.gbgout$", full.names = TRUE), method = "radix")
excl <- sub("\\.bam$", "", trimws(readLines(file.path(raw, "exclude.lst"))))
files <- files[!sub("\\.bam\\.gbgout$", "", basename(files)) %in% excl]
yc <- NULL
for (f in files) {
  t <- read.delim(f, header = FALSE, col.names = c("gene", "n"), quote = "")
  t <- t[!t$gene %in% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"), ]
  if (is.null(yc)) yc <- data.frame(gene = t$gene) else stopifnot(identical(yc$gene, t$gene))
  nm <- tolower(sub("^(WT|Snf2)_(rep[0-9]+)_.*", "\\1_\\2", basename(f)))
  yc[[nm]] <- t$n
}
ywt <- grep("^wt_", names(yc), value = TRUE)
ysn <- grep("^snf2_", names(yc), value = TRUE)
set.seed(1)
ywt <- sample(ywt); ysn <- sample(ysn) # fixture k takes the first replicates of these seeded orders

## Pasilla.
pc <- read.delim(file.path(raw, "pasilla/pasilla_gene_counts.tsv"), check.names = FALSE)
names(pc)[1] <- "gene"
pa <- read.csv(file.path(raw, "pasilla/pasilla_sample_annotation.csv"), check.names = FALSE)
pa <- data.frame(sample = sub("fb$", "", pa$file), condition = pa$condition,
                 type = ifelse(pa$type == "single-read", "single", "paired"))

## A. fumigatus.
ac <- read.delim(file.path(raw, "afum/counts.full.tsv"), check.names = FALSE)
names(ac)[1] <- "gene"
acd <- read.delim(file.path(raw, "afum/coldata.tsv"), check.names = FALSE)

write_case <- function(dir, counts, coldata, design, contrast, ref, purpose) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write.table(counts, file.path(dir, "counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(coldata, file.path(dir, "coldata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(paste(c("design", "contrast", "ref", "purpose"), c(design, contrast, ref, purpose), sep = "\t"),
             file.path(dir, "spec.tsv"))
}

# The fixture's genes: 300 at random (seed), plus up to 3 all-zero genes, up to 5 genes zero in exactly one of the two
# contrasted groups, and up to 3 genes zero in both contrasted levels but not overall; in the matrix's gene order.
subset_genes <- function(counts, coldata, contrast, seed) {
  set.seed(seed)
  m <- as.matrix(counts[, -1])
  f <- strsplit(contrast, ",")[[1]]
  g1 <- coldata$sample[coldata[[f[1]]] == f[2]]; g2 <- coldata$sample[coldata[[f[1]]] == f[3]]
  z1 <- rowSums(m[, g1, drop = FALSE]) == 0; z2 <- rowSums(m[, g2, drop = FALSE]) == 0
  all0 <- rowSums(m) == 0
  pick <- function(w, k) { w <- which(w); if (length(w) <= k) w else sample(w, k) }
  idx <- c(sample(nrow(m), 300), pick(all0, 3), pick(xor(z1, z2) & !all0, 5), pick(z1 & z2 & !all0, 3))
  counts[sort(unique(idx)), ]
}

cases <- list()
add <- function(name, counts, coldata, design, contrast, ref, purpose, seed, full = FALSE) {
  cases[[length(cases) + 1]] <<- name
  counts <- counts[, c("gene", coldata$sample)]
  write_case(file.path(outdir, "fixtures", name), subset_genes(counts, coldata, contrast, seed), coldata, design,
             contrast, ref, purpose)
  if (full) write_case(file.path(outdir, "full", name), counts, coldata, design, contrast, ref, purpose)
}
ycd <- function(wt, sn, cond = c("wt", "snf2")) {
  data.frame(sample = c(wt, sn), condition = c(rep(cond[1], length(wt)), rep(cond[2], length(sn))))
}

add("yeast_3v3", yc, ycd(ywt[1:3], ysn[1:3]), "~ condition", "condition,snf2,wt", "wt",
    "3 vs 3, m - p = 4 (trigamma prior)", 11, full = TRUE)
add("yeast_2v2", yc, ycd(ywt[1:2], ysn[1:2]), "~ condition", "condition,snf2,wt", "wt",
    "2 vs 2, m - p = 2 (simulated prior)", 12)
add("yeast_1v2", yc, ycd(ywt[1], ysn[1:2]), "~ condition", "condition,snf2,wt", "wt",
    "1 vs 2, m - p = 1 (simulated prior), a group of one", 13)
add("yeast_1v3", yc, ycd(ywt[1], ysn[1:3]), "~ condition", "condition,snf2,wt", "wt",
    "1 vs 3, m - p = 2, a group of one", 14)
add("yeast_2v5", yc, ycd(ywt[1:2], ysn[1:5]), "~ condition", "condition,snf2,wt", "wt",
    "unbalanced 2 vs 5, m - p = 5", 15)
b <- ycd(ywt[1:3], ysn[1:3]); b$batch <- c("a", "b", "a", "b", "a", "b")
add("yeast_batch_3v3", yc, b, "~ batch + condition", "condition,snf2,wt", "wt",
    "3 vs 3 plus a synthetic two-level batch, m - p = 3 (simulated prior)", 16)
x <- ycd(ywt[1:3], ysn[1:3]); x$score <- c(0.31, 1.27, 0.84, 1.62, 0.05, 0.97)
add("yeast_numeric_3v3", yc, x, "~ score + condition", "condition,snf2,wt", "wt",
    "3 vs 3 plus a synthetic numeric covariate (every sample its own Cook's cell)", 17)
l3 <- data.frame(sample = c(ywt[1:3], ywt[4:6], ysn[1:3]), group = rep(c("a", "b", "c"), each = 3))
add("yeast_3level", yc, l3, "~ group", "group,c,b", "a",
    "3-level factor (a, b: WT; c: snf2), contrast between the two non-reference levels", 18)
add("yeast_null_3v3", yc, ycd(ywt[1:3], ywt[4:6], c("g1", "g2")), "~ condition", "condition,g2,g1", "g1",
    "WT vs WT: a null comparison (expect max numRej <= 10, no filtering)", 19)
# 7 vs 7 with planted outliers (Cook's replacement and refit, phase 5): in 5 seeded genes with a mean count >= 50, one
# sample's count times 20. The full-data run below has none planted.
o <- yc[, c("gene", ywt[1:7], ysn[1:7])]
set.seed(20)
og <- sample(which(rowMeans(as.matrix(o[, -1])) >= 50), 5)
for (k in seq_along(og)) o[og[k], 1 + k] <- o[og[k], 1 + k] * 20L
add("yeast_7v7_outliers", o, ycd(ywt[1:7], ysn[1:7]), "~ condition", "condition,snf2,wt", "wt",
    paste0("7 vs 7 with planted outliers (count x20) in ", paste(o$gene[og], collapse = ","),
           ": Cook's replacement and refit"), 20)
write_case(file.path(outdir, "full", "yeast_all"), yc[, c("gene", ywt, ysn)], ycd(ywt, ysn), "~ condition",
           "condition,snf2,wt", "wt", sprintf("all good replicates, %d vs %d", length(ywt), length(ysn)))

add("pasilla_condition", pc, pa[, c("sample", "condition")], "~ condition", "condition,treated,untreated",
    "untreated", "pasilla, 3 treated vs 4 untreated", 31, full = TRUE)
add("pasilla_type_condition", pc, pa, "~ type + condition", "condition,treated,untreated", "untreated",
    "pasilla with library type (single, paired), m - p = 4", 32, full = TRUE)
add("pasilla_interaction", pc, pa, "~ type * condition", "condition,treated,untreated", "untreated",
    "interaction design: the contrast is the effect at the reference type (paired), m - p = 3", 33)

add("afum_condition", ac, acd[, c("sample", "condition")], "~ condition", "condition,Azole,GMM", "GMM",
    "A. fumigatus Azole vs GMM, kallisto counts from the page", 41, full = TRUE)
add("afum_strain_condition", ac, acd, "~ strain + condition", "condition,Azole,GMM", "GMM",
    "A. fumigatus with strain (AF293, CEA10), m - p = 4", 42, full = TRUE)

cat(sprintf("%d fixtures in %s/fixtures; full inputs in %s/full\n", length(cases), outdir, outdir))
