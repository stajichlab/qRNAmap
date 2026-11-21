# Shrinkage on the AMB counts

Input: `run_align/counts.tsv.zst` (align-mode counts, HISAT2 + featureCounts) and `samples.tsv` of the AMB run, contrast `condition,AmpB,Control`, design `~ condition`.
Commands (one per method, `METHOD` is normal, apeglm or ashr):

    qrnamap de --counts run_align/counts.tsv.zst --samples samples.tsv \
      --contrast condition,AmpB,Control --shrink METHOD --out shrink_overhang/METHOD

Versions: R 4.5.3, DESeq2 1.50.2, apeglm 1.32.0, ashr 2.2.63.
`shrink_summary.tsv` lists row counts and the sum of absolute log2 fold changes.
`shrink_correlation.tsv` lists Pearson and Spearman correlations on genes with a finite value in both columns.
The full output stays outside git.
The output directory name `shrink_overhang` has no link to the Salmon overhang runs.
