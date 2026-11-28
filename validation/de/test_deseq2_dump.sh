#!/usr/bin/env bash
# Fixture test of deseq2_dump.R, deseq2_prior_table.R and check_prior_table.R on simulated counts (needs R with DESeq2).
#   bash validation/de/test_deseq2_dump.sh        (RSCRIPT=... to choose R; default: pixi environment `de`)
# Checks: every output file and column; results() equals the stored coefficient for a two-group design; levels in C-locale
# order with the reference first; the m - p <= 3 path is taken for 2 vs 2 and the prior table reproduces its dispPriorVar.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
read -r -a R <<< "${RSCRIPT:-pixi run --manifest-path $ROOT/pixi.toml -e de Rscript}"
T="$(mktemp -d "${TMPDIR:-/tmp}/deseq2dump.XXXXXX")"; trap 'rm -rf "$T"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }
"${R[@]}" -e '
set.seed(7); a <- commandArgs(TRUE)[1]
sim <- function(n, s, lfc) { mu <- exp(runif(n, 2, 8)); disp <- 0.05 + 2 / mu
  sapply(seq_along(s), function(j) rnbinom(n, mu = mu * ifelse(s[j] == "B" & seq_len(n) <= 40, 2^lfc, 1), size = 1 / disp)) }
for (case in list(list("s33", c("A","A","A","B","B","B")), list("s22", c("B","A","B","A")))) {
  x <- sim(400, case[[2]], 2); x[5, ] <- 0
  sm <- paste0("smp", seq_along(case[[2]])); colnames(x) <- sm
  d <- file.path(a, case[[1]]); dir.create(d)
  write.table(data.frame(gene = sprintf("g%03d", 1:400), x, check.names = FALSE), file.path(d, "counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(sample = sm, grp = case[[2]]), file.path(d, "coldata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}' "$T" >/dev/null
for c in s33 s22; do
  "${R[@]}" "$HERE/deseq2_dump.R" "$T/$c/counts.tsv" "$T/$c/coldata.tsv" "~ grp" "grp,A,B" B "$T/out_$c" > "$T/log_$c" \
    || { cat "$T/log_$c"; fail "dump $c"; }
  for f in genes.tsv sizefactors.tsv hat.tsv cooks.tsv filter.tsv run.tsv messages.txt; do
    [[ -e "$T/out_$c/$f" ]] || fail "$c: no $f"
  done
done
val() { awk -F'\t' -v k="$2" '$1 == k { print $2 }' "$T/out_$1/run.tsv"; }
[[ $(val s33 levels.grp) == "B,A" ]] || fail "levels: reference first, got $(val s33 levels.grp)"
[[ $(val s33 dispPriorVar_path) == trigamma && $(val s22 dispPriorVar_path) == simulation ]] || fail "prior paths"
[[ $(val s33 m) == 6 && $(val s33 p) == 2 && $(val s33 genes) == 400 && $(val s33 genes_all_zero) == 1 ]] || fail "run.tsv counts"
"${R[@]}" -e '
g <- read.delim(file.path(commandArgs(TRUE)[1], "genes.tsv"), check.names = FALSE)
need <- c("gene", "baseMean", "allZero", "dispGeneEst", "dispGeneIter", "dispFit", "dispersion", "dispMAP", "dispOutlier",
          "Intercept", "grp_A_vs_B", "SE_grp_A_vs_B", "betaConv", "betaIter", "maxCooks", "res_log2FoldChange", "res_padj", "optimRow")
miss <- setdiff(need, names(g)); if (length(miss)) stop("missing columns: ", paste(miss, collapse = " "))
nz <- !g$allZero
stopifnot(sum(!nz) == 1, all(is.na(g$res_pvalue[!nz])))
stopifnot(isTRUE(all.equal(g$res_log2FoldChange[nz], g$grp_A_vs_B[nz], tolerance = 0)))
stopifnot(all(g$dispersion[nz] >= 1e-8 & g$dispersion[nz] <= 10))
stopifnot(sum(g$res_padj[1:40] < 0.05, na.rm = TRUE) >= 20)
h <- read.delim(file.path(commandArgs(TRUE)[1], "hat.tsv"), check.names = FALSE); stopifnot(nrow(h) == 400, ncol(h) == 7)
' "$T/out_s33" || fail "genes.tsv checks"
"${R[@]}" "$HERE/deseq2_prior_table.R" "$T/prior.tsv" >/dev/null
[[ $(wc -l < "$T/prior.tsv") -eq 601 ]] || fail "prior table rows"
"${R[@]}" "$HERE/check_prior_table.R" "$T/prior.tsv" "$T/out_s33" "$T/out_s22" | tee "$T/chk" | grep -q '^ok .*out_s22' \
  || { cat "$T/chk"; fail "the prior table does not reproduce DESeq2's dispPriorVar"; }
echo "ok: DESeq2 dump, prior table and its check"
