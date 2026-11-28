#!/usr/bin/env Rscript
# DESeq2 with default arguments on one count matrix, writing every intermediate the Rust port is checked against
# (docs/superpowers/specs/2026-11-27-deseq2-rust-design.md, "Validation" 1). R DESeq2 from the pixi environment.
#   Rscript deseq2_dump.R COUNTS.tsv COLDATA.tsv DESIGN CONTRAST REF OUTDIR
# COUNTS: `gene` column then one non-negative integer column per sample (non-integers are an error: round first).
# COLDATA: column `sample` (the count columns) plus covariates; columns R reads as numbers stay numeric.
# DESIGN: an R formula ("~ condition", "~ batch + condition", "~ type * condition").
# CONTRAST: factor,numerator,denominator (as `qrnamap de --contrast`). REF: the reference level of that factor.
# Factor levels: C-locale (byte) order, with REF first for the contrast factor; written to run.tsv so the port uses the
# same model matrix.
# OUTDIR gets (numbers with 17 significant digits, NA as NA):
#   genes.tsv        per gene: all of mcols(dds) (baseMean, baseVar, allZero, dispGeneEst, dispGeneIter, dispFit,
#                    dispersion, dispIter, dispOutlier, dispMAP, coefficients, SE_*, WaldStatistic_*, WaldPvalue_*,
#                    betaConv, betaIter, deviance, maxCooks, replace), then res_* (results() columns), then optimRow
#                    (IRLS reached maxit, so DESeq2 used its optim fallback; betaIter keeps the IRLS count)
#   sizefactors.tsv  sample, sizeFactor
#   hat.tsv, cooks.tsv  gene × sample: hat diagonals (assay H) and Cook's distances (assay cooks)
#   filter.tsv       theta, numRej, lowess fit (independent filtering)
#   run.tsv          key, value: version, design, contrast, levels, m, p, coefficients, trend type and coefficients,
#                    varLogDispEsts, dispPriorVar and its path, Cook's cutoff, filter threshold and theta, counts
#   messages.txt     DESeq2's messages and warnings
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 6) stop("usage: deseq2_dump.R COUNTS COLDATA DESIGN CONTRAST REF OUTDIR")
contrast <- strsplit(a[4], ",", fixed = TRUE)[[1]]
if (length(contrast) != 3) stop("CONTRAST must be factor,numerator,denominator")
ref <- a[5]
out <- a[6]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

raw <- read.delim(a[1], row.names = 1, check.names = FALSE)
counts <- as.matrix(raw)
if (any(is.na(counts)) || any(counts < 0) || any(counts != round(counts))) stop("counts must be non-negative integers")
storage.mode(counts) <- "integer"

cd <- read.delim(a[2], check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample" %in% names(cd)) stop("COLDATA needs a sample column")
rownames(cd) <- cd$sample
if (!setequal(cd$sample, colnames(counts))) stop("COLDATA samples differ from the count columns")
cd <- cd[colnames(counts), , drop = FALSE]
levels_out <- character()
for (v in setdiff(names(cd), "sample")) {
  if (is.numeric(cd[[v]])) next
  lv <- sort(unique(cd[[v]]), method = "radix") # radix sort is C-locale; R's default sort() depends on the locale
  if (v == contrast[1]) {
    if (!ref %in% lv) stop("REF '", ref, "' is not a level of ", v)
    lv <- c(ref, setdiff(lv, ref))
  }
  cd[[v]] <- factor(cd[[v]], levels = lv)
  levels_out[v] <- paste(lv, collapse = ",")
}

msgs <- character()
keep <- function(prefix) function(cnd) {
  msgs <<- c(msgs, paste0(prefix, sub("\n$", "", conditionMessage(cnd))))
  invokeRestart(if (inherits(cnd, "warning")) "muffleWarning" else "muffleMessage")
}
dds <- withCallingHandlers({
  dds <- DESeqDataSetFromMatrix(counts, cd, as.formula(a[3]))
  DESeq(dds)
}, message = keep(""), warning = keep("warning: "))
res <- withCallingHandlers(results(dds, contrast = contrast), message = keep(""), warning = keep("warning: "))

fmt <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE")))
  if (is.numeric(x)) return(ifelse(is.na(x), "NA", sprintf("%.17g", x)))
  ifelse(is.na(x), "NA", as.character(x))
}
put <- function(df, file) {
  df[] <- lapply(df, fmt)
  write.table(df, file.path(out, file), sep = "\t", quote = FALSE, row.names = FALSE)
}

mc <- as.data.frame(mcols(dds))
mc <- mc[, vapply(mc, function(x) is.numeric(x) || is.logical(x), TRUE), drop = FALSE]
rs <- as.data.frame(res)[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
names(rs) <- paste0("res_", names(rs))
maxit <- 100L # nbinomWaldTest default
genes <- data.frame(gene = rownames(dds), mc, rs, optimRow = !is.na(mc$betaIter) & mc$betaIter >= maxit,
                    check.names = FALSE)
put(genes, "genes.tsv")
put(data.frame(sample = colnames(dds), sizeFactor = sizeFactors(dds)), "sizefactors.tsv")
for (key in c("H", "cooks")) {
  put(data.frame(gene = rownames(dds), assays(dds)[[key]], check.names = FALSE), c(H = "hat.tsv", cooks = "cooks.tsv")[[key]])
}

md <- metadata(res)
lf <- md$lowessFit
filt <- data.frame(theta = md$filterNumRej$theta, numRej = md$filterNumRej$numRej,
                   lowess = if (is.null(lf) || length(lf$y) != nrow(md$filterNumRej)) NA else lf$y)
put(filt, "filter.tsv")

mm <- model.matrix(design(dds), colData(dds))
m <- nrow(mm); p <- ncol(mm)
df <- dispersionFunction(dds)
trend <- if (any(grepl("local regression fit", msgs))) "local" else if (any(grepl("mean of gene-wise", msgs))) "mean" else "parametric"
co <- attr(df, "coefficients")
nz <- !mcols(dds)$allZero
run <- c(
  deseq2_version = as.character(packageVersion("DESeq2")),
  design = a[3], contrast = a[4], ref = ref,
  setNames(levels_out, paste0("levels.", names(levels_out))),
  samples = paste(colnames(dds), collapse = ","),
  m = m, p = p, coefficients = paste(resultsNames(dds), collapse = ","),
  model_matrix_columns = paste(colnames(mm), collapse = ","),
  genes = nrow(dds), genes_all_zero = sum(!nz),
  trend = trend,
  trend_asympDisp = if (is.null(co)) NA else fmt(co[["asympDisp"]]),
  trend_extraPois = if (is.null(co)) NA else fmt(co[["extraPois"]]),
  varLogDispEsts = fmt(attr(df, "varLogDispEsts")),
  dispPriorVar = fmt(attr(df, "dispPriorVar")),
  dispPriorVar_path = if (m - p <= 3 && m > p) "simulation" else if (m > p) "trigamma" else "none",
  cooks_cutoff = fmt(qf(0.99, p, m - p)),
  filter_alpha = fmt(md$alpha),
  filter_threshold = fmt(md$filterThreshold),
  filter_theta = fmt(md$filterTheta),
  n_disp_outlier = sum(mcols(dds)$dispOutlier, na.rm = TRUE),
  n_gene_grid = sum(nz & (mcols(dds)$dispGeneIter %in% c(1, 100)) & mcols(dds)$dispGeneEst > 1e-7, na.rm = TRUE),
  n_optim_rows = sum(genes$optimRow),
  n_cooks_flagged = sum(nz & !is.na(res$baseMean) & is.na(res$pvalue) & !is.na(res$stat)),
  n_replaced = if (is.null(mcols(dds)$replace)) 0 else sum(mcols(dds)$replace, na.rm = TRUE),
  n_filtered = sum(nz & !is.na(res$pvalue) & is.na(res$padj)),
  n_padj_0.05 = sum(res$padj < 0.05, na.rm = TRUE)
)
writeLines(paste(names(run), run, sep = "\t"), file.path(out, "run.tsv"))
writeLines(msgs, file.path(out, "messages.txt"))
cat(sprintf("%s: m %d, p %d, trend %s, prior %s, padj < 0.05: %s\n", out, m, p, trend, run[["dispPriorVar_path"]],
            run[["n_padj_0.05"]]))
