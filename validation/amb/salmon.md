# AMB: Salmon compared with kallisto and with the authors' GSE80357 tables

Run on 2026-11-20 (US Pacific time) on the UCR HPCC cluster (SLURM 24.11.1), in the worktree `salmon`
(branch `salmon`). All numbers in this file come from the runs described here. There is no ground
truth. kallisto, Salmon and the authors' pipeline are all pipeline outputs. No criterion for "good"
agreement was set, so this file lists numbers only.

Data, samples and strandedness are the same as in `validation/amb/README.md` (6 paired-end samples,
8,405,534 to 10,302,142 pairs; strand-check result `unstranded`). The GEO tables and their mapping to
our samples are from `validation/geo/README.md`. `$DATA` below is
`/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb`.

## Files

- `run_amb_salmon.sh`: the driver. Stages `gen`, `submit`, `compare` (see the script header).
- `results_salmon/mode_agreement.tsv`, `results_salmon/de_overlap.tsv`: `compare_modes.R` with Salmon in
  the QUANT slot and kallisto in the ALIGN slot. In these two files, "quant" is Salmon and "align" is
  kallisto.
- `results_salmon/geo/*.tsv`: `validation/geo/compare_geo.R` with Salmon in the QUANT slot and kallisto in
  the ALIGN slot. In these files, "quant" is Salmon and "align" is kallisto. No row in them is HISAT2.
- `results_salmon/mapping_summary.tsv`: the decompressed `mapping_summary.tsv.zst` of the Salmon run.
- `results_salmon/jobs.tsv`: `sacct` output of the four jobs of the run and the compare step.
- `collect_salmon.py` and `results_salmon/per_sample.tsv`: per-sample values from `meta_info.json`
  (`salmon_version`, `num_processed`, `num_mapped`, `percent_mapped`) and the per-sample seconds.
- `zero_genes_by_length.py` and `results_salmon/zero_genes_by_length.tsv`, `short_transcripts.tsv`,
  `zero_genes_examples.tsv`: the genes that are 0 in Salmon and nonzero in kallisto, by transcript length.
- `salmon_subset.sh` and `results_salmon/subset_Control_1.tsv`: Salmon option variants and kallisto on
  the first 500,000 pairs of Control_1 (job 32329651).

Outputs outside the repository: `$DATA/run_quant_salmon/` (run directory), `$DATA/geo/de_salmon_rerun/`,
logs `$DATA/logs/gen_salmon.log`, `$DATA/logs/submit_salmon.log`,
`$DATA/logs/salmon_compare-32328937.out`, `$DATA/run_quant_salmon/slurm/logs/`.

## Commands

```bash
W=/bigdata/stajichlab/jstajich/projects/qRNAmap/.claude/worktrees/salmon
D=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb
export PATH=$W/.pixi/envs/default/bin:$PATH
CARGO_TARGET_DIR=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-salmon \
  cargo build --release -p qrnamap-cli                 # commit ed9d5b6
REPO=$W bash $W/validation/amb/run_amb_salmon.sh gen      # login node
REPO=$W bash $W/validation/amb/run_amb_salmon.sh submit   # login node, fresh login shell
env -i HOME="$HOME" USER="$USER" PATH=/usr/bin:/bin bash -l -c \
  "sbatch --parsable -p short -c 1 --mem=4G --time=00:30:00 -J salmon_compare \
   --output=$D/logs/salmon_compare-%j.out --wrap 'REPO=$W bash $W/validation/amb/run_amb_salmon.sh compare'"
# per-sample table and zero-count genes (login node, read only on $D)
python -I $W/validation/amb/collect_salmon.py --run-dir $D/run_quant_salmon --sheet $D/samples.tsv \
  --task-start 2026-11-20T17:13:47 --out $W/validation/amb/results_salmon/per_sample.tsv
python -I $W/validation/amb/zero_genes_by_length.py --sheet $D/samples.tsv --salmon-dir $D/run_quant_salmon \
  --kallisto-dir $D/run_quant --authors $D/geo/prep/authors_counts_amb.tsv --out-dir $W/validation/amb/results_salmon
# subset diagnostic (job 32329651, node c23, 2 CPUs, 0:57, MaxRSS 460,840 KB)
env -i HOME="$HOME" USER="$USER" PATH=/usr/bin:/bin bash -l -c \
  "sbatch --parsable -p short -c 2 --mem=8G --time=00:30:00 -J salmon_subset \
   --output=$D/logs/salmon_subset-%j.out --wrap 'REPO=$W bash $W/validation/amb/salmon_subset.sh Control_1 500000'"
```

When the run and compare jobs ran, `run_amb_salmon.sh` gave `REPO` a default equal to `$W`. The
script now requires `REPO`, so the commands above set it.

The `gen` stage runs:
`qrnamap run --samples $DATA/samples.tsv --reference GCF_000146045.2_R64_genomic.fna.gz --gff
GCF_000146045.2_R64_genomic.gff.gz --mode quant --quantifier salmon --out $DATA/run_quant_salmon
--strandedness unstranded --threads 1 --contrast condition,AmpB,Control --slurm --samples-per-task 6
--partition short --time 02:00:00`. It uses the default `--mem-gb 8` and no `--tmpdir`, so the jobs put
temporary files in node-local `$SCRATCH`. The kallisto run (`validation/amb/README.md`, E2) used the
same options without `--quantifier`, and gave `--mem-gb 8` explicitly.

Versions:

- qrnamap: built from commit `ed9d5b6`.
- Salmon: 2.8.0. Every `aux_info/meta_info.json` gives this `salmon_version` (`results_salmon/per_sample.tsv`).
- R 4.5.3 and DESeq2 1.50.2, from `deseq2.versions.txt`.
- DESeq2 settings: design `~ condition`, contrast `condition,AmpB,Control`, Wald test, no shrinkage.
- kallisto side: the earlier run `$DATA/run_quant`. It used kallisto 0.52.0 and a binary from commit
  `ff611ae` (`validation/amb/README.md`). I did not run kallisto again.

## Jobs

Elapsed and MaxRSS are from `sacct` (the `.batch` step). All jobs ended `COMPLETED 0:0`. Every
qrnamap job asked for 1 CPU and `--mem=8G`.

| job | step | node | AllocCPUS | elapsed | MaxRSS |
|---|---|---|---|---|---|
| 32328232 | prepare (transcripts, `salmon index`) | r43 | 2 | 1:58 | 449,860 KB |
| 32328233_0 | sample array, task 0 (6 samples) | x01 | 2 | 17:25 | 84,344 KB |
| 32328234 | final (matrix, DESeq2) | c19 | 1 | 0:50 | 835,568 KB |
| 32328937 | compare (`--mem=4G`) | r43 | 2 | 0:51 | 827,516 KB |

For reference, the kallisto run (`validation/amb/README.md`, E2): prepare 1:44 (182,448 KB, r43),
sample task 23:29 (118,460 KB, r41), final 0:29 (870,664 KB, x05). The two sample tasks ran on
different nodes (r41 and x01) at different times. This is not a controlled speed comparison.

Size on disk, read from files outside the repository with `du -sb`: `$DATA/run_quant_salmon/salmon_index/`
216,730,680 bytes (211 MB with `du -sh`); `$DATA/run_quant/transcripts.kidx` 5,835,107 bytes.

## Per sample

Source: `results_salmon/per_sample.tsv` (`collect_salmon.py`). Pairs, `num_mapped` and mapped % come
from `meta_info.json`. Seconds are derived from the mtimes of `SAMPLE.genecounts.tsv.zst`, as in
`validation/amb/README.md`: the difference to the previous sample, or to the task's sacct Start for the
first sample. They are not timer values. The kallisto pseudoaligned % is from `validation/amb/README.md`.

| sample | pairs | Salmon mapped % | Salmon `num_mapped` | Salmon s | Salmon s per million pairs | kallisto pseudoaligned % |
|---|---|---|---|---|---|---|
| Control_1 | 9,625,471 | 87.54 | 8,425,674 | 192 | 19.9 | 93.6 |
| Control_2 | 10,302,142 | 90.45 | 9,318,680 | 201 | 19.5 | 94.6 |
| Control_3 | 9,997,148 | 90.26 | 9,023,103 | 171 | 17.1 | 94.7 |
| AmpB_1 | 8,405,534 | 89.34 | 7,509,631 | 136 | 16.1 | 94.0 |
| AmpB_2 | 9,755,935 | 88.04 | 8,588,953 | 155 | 15.9 | 93.6 |
| AmpB_3 | 10,175,785 | 89.22 | 9,079,039 | 191 | 18.8 | 94.1 |

- Salmon `num_processed` equals the ENA pair count for all 6 samples (the qrnamap check passed).
- Whole sample task: 1,045 s for 58,262,015 pairs = 17.9 s per million pairs (kallisto task: 24.2).
- No line starting with `warning:` appeared in the SLURM logs. Read from files outside the repository:
  `grep -c '^warning:' $DATA/run_quant_salmon/slurm/logs/*.out` gave 0 for each of the 3 logs.

## (a) Salmon compared with kallisto

`compare_modes.R $DATA/run_quant_salmon $DATA/run_quant`. Both matrices have 6,477 gene rows. Row and
column names are identical and in the same order. Spearman: R `cor(method = "spearman")`. Pooled =
6 samples x 6,477 genes = 38,862 values. `results_salmon/mode_agreement.tsv`:

| sample | Spearman, all genes | Spearman, genes nonzero in either (n) | Pearson of log2(count+1) | total Salmon | total kallisto | zero genes Salmon | zero genes kallisto |
|---|---|---|---|---|---|---|---|
| Control_1 | 0.994519 | 0.994444 (6,371) | 0.978292 | 8,425,634 | 9,005,877 | 411 | 124 |
| Control_2 | 0.995464 | 0.995373 (6,349) | 0.984154 | 9,318,622 | 9,744,045 | 416 | 136 |
| Control_3 | 0.995372 | 0.995293 (6,339) | 0.983119 | 9,023,021 | 9,462,768 | 430 | 152 |
| AmpB_1 | 0.994639 | 0.994523 (6,343) | 0.982641 | 7,509,522 | 7,905,110 | 432 | 144 |
| AmpB_2 | 0.994674 | 0.994579 (6,362) | 0.981232 | 8,588,901 | 9,134,449 | 420 | 126 |
| AmpB_3 | 0.994712 | 0.994627 (6,370) | 0.979950 | 9,078,969 | 9,580,410 | 413 | 119 |
| pooled | 0.994919 | 0.994830 (38,134) | 0.981558 | 51,944,669 | 54,832,659 | 2,522 | 801 |

Salmon has more zero-count genes than kallisto. Per sample, the ratio Salmon / kallisto is 2.8 to 3.5
(lowest Control_3, 430/152 = 2.83; highest AmpB_3, 413/119 = 3.47). Pooled, it is 2,522 against 801,
that is 3.1 times.

### Zero in Salmon, nonzero in kallisto: transcript length (measured)

Source: `results_salmon/zero_genes_by_length.tsv`, `short_transcripts.tsv`, `zero_genes_examples.tsv`
(`zero_genes_by_length.py`). The read length is 101 nt in all 6 samples (maximum over the first 10,000
records of each mate file). A gene's length class is the length of its longest transcript.

| sample | genes 0 in Salmon, >0 in kallisto | longest transcript <= 101 nt | 102-150 | 151-300 | >300 | kallisto reads in these genes |
|---|---|---|---|---|---|---|
| Control_1 | 305 | 278 | 7 | 3 | 17 | 41,353 |
| Control_2 | 288 | 267 | 7 | 2 | 12 | 31,611 |
| Control_3 | 292 | 263 | 6 | 5 | 18 | 30,113 |
| AmpB_1 | 298 | 266 | 7 | 4 | 21 | 17,369 |
| AmpB_2 | 305 | 279 | 3 | 5 | 18 | 27,411 |
| AmpB_3 | 306 | 282 | 7 | 4 | 13 | 26,814 |

- 263 to 282 of these genes per sample (89 to 93%) have no transcript longer than the reads.
- Of the 6,478 transcripts, 353 are 101 nt or shorter. Salmon gives `NumReads` 0 to 316 to 324 of them
  per sample. Of the 411 transcripts shorter than 150 nt, 326 to 333 have `NumReads` 0.
- RPL41A (YDL184C) and RPL41B (YDL133C-A) each have one 78-nt transcript. Salmon: 0 reads in all 6
  samples. kallisto: 7,690 to 19,722 (RPL41A) and 7,877 to 19,470 (RPL41B). Authors' table: 9,941 to
  23,651 and 9,696 to 23,266. In Control_1: kallisto 19,722 and 19,470; authors 23,651 and 23,266.
- YOR302W (78 nt): Salmon 0 in all 6 samples; kallisto 127 to 635; authors 92 to 512.

Hypothesis (not proven): under its default selective-alignment scoring, Salmon does not assign reads
that overhang the end of a transcript shorter than the read. A 101-nt read cannot fit inside a 78-nt
transcript. The subset run below supports this hypothesis.

Subset run (`results_salmon/subset_Control_1.tsv`, job 32329651): the first 500,000 pairs of Control_1,
one sample only, the index of `$DATA/run_quant_salmon`, `-l IU`, `-p 2`. NumReads of RPL41A + RPL41B:

| variant | RPL41A + RPL41B | YOR302W | mapped % |
|---|---|---|---|
| Salmon defaults (as qrnamap runs it) | 0 | 0 | 87.9636 |
| Salmon `--useEM` | 0 | 0 | 87.9636 |
| Salmon `--softclipOverhangs` | 163.0 | 1 | 88.1294 |
| Salmon `--softclip` | 169.0 | 1 | 88.9424 |
| Salmon `--softclipOverhangs --minScoreFraction 0.5` | 742.0 | 1 | 89.7554 |
| kallisto 0.52.0 (unstranded) | 2,101.0 | 8 | 93.8 (pseudoaligned) |

The split between RPL41A and RPL41B differs between the variants; the table gives the sum. No Salmon
variant reached the kallisto count on this subset. qrnamap passes none of these tuning flags. Changing
the Salmon flags is an open user decision (see README, "Not in this version"); the Salmon spec lists flag
tuning as a non-goal. I did not test `--noLengthCorrection` or other flags.

DESeq2, `padj < 0.05`, NA padj excluded (`results_salmon/de_overlap.tsv`):

| measure | value |
|---|---|
| genes in each result | 6,477 (Salmon), 6,477 (kallisto) |
| padj NA | 370 (Salmon), 335 (kallisto) |
| padj < 0.05, Salmon | 2,005 |
| padj < 0.05, kallisto | 2,042 |
| in both | 1,962 |
| Salmon only | 43 |
| kallisto only | 80 |
| Jaccard (both / either) | 0.941007 |
| in both with the same log2FC sign | 1,962 of 1,962 |
| Spearman of log2FC, genes with a log2FC in both | 0.974197 |

For comparison, kallisto vs HISAT2 + featureCounts on the same data (`validation/amb/README.md`):
pooled Spearman 0.964905, DE Jaccard 0.929216.

## (b) Salmon and kallisto compared with the authors' tables

`compare_geo.R` (unchanged) with Salmon in the QUANT slot and kallisto in the ALIGN slot. The
authors' count table is `$DATA/geo/prep/authors_counts_amb.tsv` from the GEO comparison. Common genes
(in our matrix and in the authors' counts): 6,008. Definitions are in `validation/geo/README.md`.

Check: two kallisto rows of this run were compared field by field with the earlier GEO results
(`validation/geo/results/`). The pooled `authors(x)_vs_align(y)` count row equals the pooled
`authors(x)_vs_quant(y)` row there. The `authors_published(x)_vs_deseq2_align(y)` DE row equals the
`authors_published(x)_vs_deseq2_quant(y)` row there. I did not compare the other rows.

Check of the DESeq2 step: `qrnamap de` from `ed9d5b6` on the Salmon counts, against the
`deseq2.results.csv` of the SLURM final job: 6,477 genes; max |log2FC difference| 2.6e-10; NA padj in the
same genes; max |padj difference| 5.1e-08. Read from a file outside the repository:
`$DATA/logs/salmon_compare-32328937.out`, line 131 ("rerun check").

### Counts (`results_salmon/geo/counts_agreement.tsv`, pooled rows, 36,048 values)

The ratio columns are (ours + 1) / (authors + 1) over genes with an authors' count of at least 10
(34,094 values).

| comparison | Spearman | Pearson log2(x+1) | total authors | total ours | total ratio ours/authors | median ratio | IQR of ratio |
|---|---|---|---|---|---|---|---|
| authors vs Salmon | 0.941675 | 0.876650 | 46,263,035 | 51,129,364 | 1.105188 | 1.016304 | 0.069083 |
| authors vs kallisto | 0.946995 | 0.890594 | 46,263,035 | 53,986,917 | 1.166956 | 1.058282 | 0.036250 |
| kallisto (x) vs Salmon (y), 6,008 genes | 0.994663 | 0.981354 | 53,986,917 | 51,129,364 | 0.947070 | 0.957198 | 0.045904 |

Per sample, authors vs Salmon Spearman: 0.939495 to 0.943530. Authors vs kallisto: 0.944480 to
0.949078.

### DE (`results_salmon/geo/de_agreement.tsv`)

Significant: FDR or padj < 0.05. The authors' DE table has 5,456 genes; 5,407 of them are in our
results.

| comparison | common genes | Spearman log2FC | Pearson log2FC | sig authors (common) | sig ours (common) | both | authors only | ours only | Jaccard | both, same sign |
|---|---|---|---|---|---|---|---|---|---|---|
| authors published vs DESeq2 Salmon | 5,407 | 0.986613 | 0.968686 | 1,816 | 1,922 | 1,638 | 178 | 284 | 0.780000 | 1,638 |
| authors published vs DESeq2 kallisto | 5,407 | 0.988440 | 0.956567 | 1,816 | 1,955 | 1,657 | 159 | 298 | 0.783822 | 1,657 |
| DESeq2 on authors' counts vs DESeq2 Salmon | 6,008 | 0.954430 | 0.819001 | 1,941 | 1,980 | 1,846 | 95 | 134 | 0.889639 | 1,846 |
| DESeq2 on authors' counts vs DESeq2 kallisto | 6,008 | 0.962264 | 0.889690 | 1,941 | 2,013 | 1,874 | 67 | 139 | 0.900962 | 1,874 |
| DESeq2 on authors' counts vs DESeq2 Salmon, authors' DE genes only | 5,407 | 0.991166 | 0.968418 | 1,919 | 1,922 | 1,833 | 86 | 89 | 0.912849 | 1,833 |
| DESeq2 on authors' counts vs DESeq2 kallisto, authors' DE genes only | 5,407 | 0.993407 | 0.957969 | 1,919 | 1,955 | 1,859 | 60 | 96 | 0.922581 | 1,859 |

The last two rows use only the 5,407 genes of the authors' DE table that are in our results. In the
"DESeq2 on authors' counts" rows, "authors" means our DESeq2 step run on the authors' count
table (`$DATA/geo/de_authors/`, from the GEO comparison).

Sign check (`results_salmon/geo/sign_convention.tsv`): Salmon DESeq2 log2FC vs the authors'
`AMB_only_logFC`: 5,404 genes, Spearman 0.986613; among the 1,814 genes with authors' FDR < 0.05 and a
Salmon log2FC, 1,814 have the same sign.

## Summary of the measured differences

- Salmon mapped 87.54 to 90.45% of pairs; kallisto pseudoaligned 93.6 to 94.7% on the same samples.
- Gene counts, Salmon vs kallisto: pooled Spearman 0.994919; Salmon totals are 0.947 of kallisto totals
  (pooled over the 6,477 genes: 51,944,669 / 54,832,659).
- Zero-count gene values: Salmon 2,522 and kallisto 801, pooled over the 6 samples. Per sample, 263 to
  282 of the genes that are 0 in Salmon and nonzero in kallisto have no transcript longer than the
  101-nt reads. RPL41A and RPL41B (78 nt) get 0 reads in Salmon and 7,690 to 19,722 in kallisto.
- DE, Salmon vs kallisto: Jaccard 0.941007 (1,962 genes in both, 43 Salmon only, 80 kallisto only).
- Against the authors' counts: Spearman 0.941675 (Salmon) and 0.946995 (kallisto); median count ratio
  ours/authors 1.016 (Salmon) and 1.058 (kallisto).
- Against the authors' published DE: Jaccard 0.780000 (Salmon) and 0.783822 (kallisto).
- Sample task time: 17.9 s per million pairs (Salmon, node x01) and 24.2 (kallisto, node r41). Different
  nodes and times.

## Not measured

- Salmon with genome decoys, bootstrap or Gibbs sampling, or bias correction (not supported).
- Salmon with more than 1 thread inside qrnamap (the subset diagnostic used `-p 2`).
- Salmon on the yeast PRJNA674653 data of `validation/README.md`.
- Whether Salmon flags such as `--softclipOverhangs` or `--minScoreFraction` change the whole-run results.
  Only the 500,000-pair subset of one sample was run with them.
- Salmon memory at `--mem-gb` below 8; the sample task MaxRSS was 84,344 KB.
- Run time of Salmon and kallisto on the same node.
