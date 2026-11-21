# deseq2_ref.R COUNTS_TSV OUT_PREFIX : DESeq2 Wald test, design ~ condition, contrast AmpB vs Control.
# Same call sequence as the generated script of crates/qrnamap-core/src/de.rs (no shrinkage).
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(DESeq2))
counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
storage.mode(counts) <- "integer"
coldata <- data.frame(condition = factor(sub("_[0-9]+$", "", colnames(counts)), levels = c("Control", "AmpB")),
                      row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds, test = "Wald", quiet = TRUE)
res <- results(dds, contrast = c("condition", "AmpB", "Control"))
out <- data.frame(gene = rownames(res), as.data.frame(res), check.names = FALSE)
write.csv(out, paste0(args[2], ".results.csv"), row.names = FALSE)
writeLines(c(paste("DESeq2", as.character(packageVersion("DESeq2"))), paste("R", R.version.string),
             "design: ~ condition", "contrast: condition,AmpB,Control",
             paste("genes", nrow(counts), "samples", ncol(counts))), paste0(args[2], ".versions.txt"))
