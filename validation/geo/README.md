# AMB: qRNAmap compared with the authors' GSE80357 tables

Run on 2026-11-20 (US Pacific time) on the UCR HPCC cluster (SLURM 24.11.1), in the worktree
`geo-compare` (branch `geo-compare`, created at `d64f917`). No CLI code was changed. All numbers in
this file come from the runs described here. There is no criterion for "good" or "bad" agreement,
and no ground truth: both sides are pipeline outputs. The tables list numbers and caveats only.

This builds on `validation/amb/README.md` (data, samples, run directories, strand check result
`unstranded`). `$DATA` below is `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb`.

## Files

- `fetch_geo.sh GEO_DIR`: downloads the three series tables, the 18 GEO sample records (text) and
  the ENA SRX-to-SRR records; writes `GEO_DIR/download/provenance.tsv` (size, md5, dates).
- `prepare_geo.py` (run with the pixi `python -I`): maps the authors' columns to our samples, writes
  the authors' AMB count matrix with our sample names, and the gene id overlap tables.
- `geo_job.sh REPO DATA QRNAMAP`: the SLURM job body. Runs `qrnamap de` on the authors' counts, runs
  `qrnamap de` again on our quant counts (a check of the binary), then `compare_geo.R`.
- `compare_geo.R`: counts agreement, DE agreement, sign-convention checks, checks of the authors'
  tables.
- `results/`: `download_provenance.tsv`, `srx_to_srr.tsv`, `column_mapping.tsv`,
  `gene_id_overlap.tsv`, `unmatched_summary.tsv`, `unmatched_top.tsv`, `counts_agreement.tsv`,
  `de_agreement.tsv`, `sign_convention.tsv`, `sign_convention_top10.tsv`,
  `authors_table_checks.tsv`, `compare_sessioninfo.txt`.

Outputs outside the repository: `$DATA/geo/download/` (the three files), `$DATA/geo/gsm/`
(GEO sample records), `$DATA/geo/gse80357_series_full.txt`, `$DATA/geo/prep/`
(`authors_counts_amb.tsv`, `unmatched_ours.tsv`, `unmatched_theirs.tsv`), `$DATA/geo/de_authors/`
(DESeq2 on the authors' counts), `$DATA/geo/de_quant_rerun/`, logs `$DATA/logs/geo_compare-JOBID.out`.

## Commands

```bash
W=/bigdata/stajichlab/jstajich/projects/qRNAmap/.claude/worktrees/geo-compare
D=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb
Q=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-harden/release/qrnamap
bash $W/validation/geo/fetch_geo.sh $D/geo                                  # login node
PATH=$W/.pixi/envs/default/bin:$PATH python -I $W/validation/geo/prepare_geo.py \
  --geo-dir $D/geo --srx-srr $D/geo/srx_to_srr.tsv --sheet $D/samples.tsv \
  --quant-counts $D/run_quant/counts.tsv.zst --align-counts $D/run_align/counts.tsv.zst \
  --out-dir $D/geo/prep --results-dir $W/validation/geo/results              # login node
env -i HOME="$HOME" USER="$USER" PATH=/usr/bin:/bin bash -l -c \
  "sbatch --parsable -p short -c 1 --mem=4G --time=00:30:00 -J geo_compare \
   --output=$D/logs/geo_compare-%j.out --wrap 'bash $W/validation/geo/geo_job.sh $W $D $Q'"
```

The downloads were first run by hand with the same URLs (15:16, 2026-11-20); `fetch_geo.sh` was
then run once (15:27) and gave the same md5 sums. The job was submitted from a login shell so that it
has `$SCRATCH` (see `validation/amb/README.md`, D3).

Jobs (all `short`, 1 CPU requested, `--mem=4G`, node c17, `COMPLETED 0:0`):

| job | elapsed | MaxRSS | note |
|---|---|---|---|
| 32322745 | 1:14 | 871,628 KB | first run; the check table had NA values (93 NA rows in the LF_only columns, see below) |
| 32322830 | 1:03 | 859,668 KB | checks with NA rows removed |
| 32322903 | 1:13 | 877,212 KB | final: adds the DE rows restricted to the authors' DE genes; all tables in `results/` are from this job |

The printed `counts_agreement.tsv` section of the three logs is byte-identical (same md5 of the
section). The printed `de_agreement.tsv` was extended in the last job (its column layout changed), so
I did not compare it byte for byte across the jobs.

Versions: `qrnamap` built from `d64f917` (`cargo build --release -p qrnamap-cli`,
`CARGO_TARGET_DIR=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-harden`). R 4.5.3,
DESeq2 1.50.2 (`$DATA/geo/de_authors/deseq2.versions.txt`; design `~ condition`, contrast
`condition,AmpB,Control`, Wald, no shrinkage). Python 3.15.2, zstd 1.5.7 (pixi env `default` of
this worktree). Our quant and align DESeq2 results are the files of the earlier AMB runs (binary
`ff611ae`). Check: `qrnamap de` from `d64f917` on the quant counts gave the same results (6,477
genes; max |log2FC difference| 1.0e-14; max |padj difference| 2.9e-13; NA padj in the same genes).

## 1. Data, provenance and mapping

Series page text (`acc.cgi?acc=GSE80357&targ=self&form=text&view=brief`) lists three files under
`ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80357/suppl/`. Downloaded over HTTPS from
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80357/suppl/<file>` on 2026-11-20
(`results/download_provenance.tsv`):

| file | bytes | md5 | Last-Modified | lines (with header) | columns |
|---|---|---|---|---|---|
| GSE80357_counts_of_mapped_reads.csv.gz | 236,984 | 18f5c27d74400d4dceb70c56e7d9935f | 16 Apr 2016 15:41:14 GMT | 6,693 | `oln_id` + 18 sample codes |
| GSE80357_log_fold_change_and_FDR_values.csv.gz | 294,334 | 373e129cd344fd46afe69f446d89766e | same | 5,457 | `genes` + 12 (logFC, LR, PValue, FDR for AMB_and_LF, AMB_only, LF_only) |
| GSE80357_average_cpm_per_treatment_group.csv.gz | 78,915 | 0cc7fbd7e7ced5ace9472a13d0b88fe0 | same | 5,457 | `systematic_name` + 6 groups |

All three pass `gzip -t`. CSV with CRLF line ends. The counts file has 6,692 genes, no duplicate id,
integer values. The DE file has 5,456 genes; the AMB_only columns have no NA; the four LF_only
columns have the text `NA` in 93 rows (also 93 NA rows in the CPM file columns 6 and 7).

Column mapping. The counts header has codes `1B, 2F, 2I, 3F, 3G, 4E, 4I, 5C, 5G, 8A, 9D, 9E,
10B, 10D, 10F, 10G, 10H, 10I`. Each GEO sample record (GSM2125072 to GSM2125089, text form, fetched
2026-11-20) has a first `!Sample_description` line with one of these codes in lower case (for
example GSM2125072: `!Sample_title = AMB only rep 1`, `!Sample_description = 2i`). Every one of
the 18 header codes matched exactly one GSM (upper-cased). The GSM `!Sample_relation = SRA:` line
gives the SRX; the ENA `filereport` of the SRX gives the run; the run is the FASTQ name in
`$DATA/samples.tsv`. The ENA sample titles agree with the GSM titles
(`results/srx_to_srr.tsv`, `results/column_mapping.tsv`):

| authors' column | GSM | GSM title | SRX | run | our sample |
|---|---|---|---|---|---|
| 3F | GSM2125075 | Matched control for AMB only rep 1 | SRX1709856 | SRR3396385 | Control_1 |
| 5G | GSM2125076 | Matched control for AMB only rep 2 | SRX1709857 | SRR3396386 | Control_2 |
| 8A | GSM2125077 | Matched control for AMB only rep 3 | SRX1709858 | SRR3396387 | Control_3 |
| 2I | GSM2125072 | AMB only rep 1 | SRX1709853 | SRR3396381 | AmpB_1 |
| 4E | GSM2125073 | AMB only rep 2 | SRX1709854 | SRR3396382 | AmpB_2 |
| 9E | GSM2125074 | AMB only rep 3 | SRX1709855 | SRR3396384 | AmpB_3 |

The other 12 columns (AMB + LF, LF only and their controls) were not used.

Authors' method (GEO `!Sample_data_processing`, same two lines in all 18 sample records): HTSeq
0.5.3p9, `-s no -m union -t CDS -i gene_id`, "reads found in ambiguous overlapping CDS were not
counted", SGD GFF converted to GTF, R64-1-1. The records and the series text do not name a DE
method (a search for "edger", "limma" and "deseq" in the 18 records found nothing). The DE column
names are `logFC`, `LR`, `PValue`, `FDR`. I did not determine the method from them.

### Gene ids

Ours: the 6,477 RefSeq GFF3 gene labels of the AMB runs. Theirs: the 6,692 ids of the counts file.
Exact string match (`results/gene_id_overlap.tsv`):

| set | genes | reads, sum of 6 samples |
|---|---|---|
| intersection | 6,008 | authors 46,263,035; quant 53,986,917; align 52,059,153 |
| only ours | 469 | quant 845,742; align 194,559 |
| only theirs | 684 | authors 173,320 |
| all ours | 6,477 | quant 54,832,659; align 52,253,712 |
| all theirs | 6,692 | authors 46,436,355 |

Unmatched by name class (`results/unmatched_summary.tsv`; classes are regular expressions on the id):

| side | class | genes | genes with reads (quant / align / authors) | reads (quant / align / authors) |
|---|---|---|---|---|
| only ours | SGD ncRNA name `YNC?####W/C` | 410 | 400 / 359 | 701,335 / 53,398 |
| only ours | systematic ORF name `Y??###W/C` | 31 | 28 / 25 | 142,935 / 140,233 |
| only ours | tRNA `t?(NNN)?` | 24 | 0 / 0 | 0 / 0 |
| only ours | mitochondrial `Q####` | 3 | 2 / 0 | 63 / 0 |
| only ours | other (`YNCI0012`) | 1 | 1 / 1 | 1,409 / 928 |
| only theirs | systematic ORF name | 675 | 498 | 173,320 |
| only theirs | mitochondrial `Q####` | 9 | 0 | 0 |

Largest unmatched genes (`results/unmatched_top.tsv`): ours YNCL0012C and YNCL0021C (quant 175,157
each, align 0), YDR134C (quant 118,171, align 120,439), YNCL0025C and YNCL0016C (quant 76,726 each,
align 0); theirs YDR133C (92,566), YLR154W-A (12,869), YLR154W-E (7,217), YLR154W-B (4,159),
YGR182C (2,322). The authors' 684 ids include 342 with at least 10 reads summed over the 6 samples.
I did not check whether pairs such as YDR134C (ours) and YDR133C (theirs) are the same locus.

## 2. Counts agreement

Gene set: the 6,008 common genes. Per sample. Ratio columns: median and quartiles of
(y+1)/(x+1) over the genes with at least 10 reads in the authors' table for that sample (the same
gene set in all three comparisons, also for align vs quant). Spearman: R `cor(method="spearman")`,
ties averaged. The Spearman of log2(count+1) equals the Spearman of raw counts in every row (a
monotone transform does not change ranks), so it is one column. `results/counts_agreement.tsv`.

| x vs y | sample | Spearman | Pearson log2(c+1) | total x | total y | total y/x | genes (x>=10 authors) | median (y+1)/(x+1) | Q1 | Q3 | IQR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| authors vs quant | Control_1 | 0.9456 | 0.8888 | 7,616,528 | 8,927,944 | 1.1722 | 5,681 | 1.0600 | 1.0423 | 1.0782 | 0.0359 |
| authors vs quant | Control_2 | 0.9449 | 0.8917 | 8,266,574 | 9,657,565 | 1.1683 | 5,687 | 1.0545 | 1.0398 | 1.0706 | 0.0308 |
| authors vs quant | Control_3 | 0.9445 | 0.8905 | 7,809,244 | 9,100,477 | 1.1653 | 5,682 | 1.0552 | 1.0397 | 1.0716 | 0.0319 |
| authors vs quant | AmpB_1 | 0.9489 | 0.8970 | 6,771,213 | 7,822,373 | 1.1552 | 5,669 | 1.0503 | 1.0329 | 1.0695 | 0.0366 |
| authors vs quant | AmpB_2 | 0.9486 | 0.8865 | 7,677,305 | 9,014,616 | 1.1742 | 5,687 | 1.0692 | 1.0515 | 1.0876 | 0.0362 |
| authors vs quant | AmpB_3 | 0.9491 | 0.8884 | 8,122,171 | 9,463,942 | 1.1652 | 5,688 | 1.0615 | 1.0444 | 1.0788 | 0.0344 |
| authors vs quant | pooled | 0.9470 | 0.8906 | 46,263,035 | 53,986,917 | 1.1670 | 34,094 | 1.0583 | 1.0410 | 1.0773 | 0.0362 |
| authors vs align | Control_1 | 0.9762 | 0.9594 | 7,616,528 | 8,568,674 | 1.1250 | 5,681 | 1.0586 | 1.0478 | 1.0746 | 0.0269 |
| authors vs align | Control_2 | 0.9761 | 0.9589 | 8,266,574 | 9,335,281 | 1.1293 | 5,687 | 1.0633 | 1.0536 | 1.0775 | 0.0239 |
| authors vs align | Control_3 | 0.9746 | 0.9563 | 7,809,244 | 8,771,984 | 1.1233 | 5,682 | 1.0594 | 1.0494 | 1.0748 | 0.0254 |
| authors vs align | AmpB_1 | 0.9772 | 0.9588 | 6,771,213 | 7,575,850 | 1.1188 | 5,669 | 1.0598 | 1.0483 | 1.0762 | 0.0279 |
| authors vs align | AmpB_2 | 0.9763 | 0.9583 | 7,677,305 | 8,674,992 | 1.1300 | 5,687 | 1.0690 | 1.0570 | 1.0855 | 0.0285 |
| authors vs align | AmpB_3 | 0.9767 | 0.9578 | 8,122,171 | 9,132,372 | 1.1244 | 5,688 | 1.0652 | 1.0551 | 1.0803 | 0.0252 |
| authors vs align | pooled | 0.9762 | 0.9583 | 46,263,035 | 52,059,153 | 1.1253 | 34,094 | 1.0627 | 1.0514 | 1.0787 | 0.0272 |
| align vs quant (reference) | Control_1 | 0.9673 | 0.9227 | 8,568,674 | 8,927,944 | 1.0419 | 5,681 | 1.0026 | 0.9890 | 1.0118 | 0.0229 |
| align vs quant (reference) | Control_2 | 0.9674 | 0.9309 | 9,335,281 | 9,657,565 | 1.0345 | 5,687 | 0.9935 | 0.9807 | 1.0000 | 0.0193 |
| align vs quant (reference) | Control_3 | 0.9671 | 0.9320 | 8,771,984 | 9,100,477 | 1.0374 | 5,682 | 0.9978 | 0.9846 | 1.0048 | 0.0202 |
| align vs quant (reference) | AmpB_1 | 0.9699 | 0.9351 | 7,575,850 | 7,822,373 | 1.0325 | 5,669 | 0.9931 | 0.9795 | 1.0007 | 0.0212 |
| align vs quant (reference) | AmpB_2 | 0.9697 | 0.9270 | 8,674,992 | 9,014,616 | 1.0391 | 5,687 | 1.0016 | 0.9877 | 1.0105 | 0.0229 |
| align vs quant (reference) | AmpB_3 | 0.9701 | 0.9279 | 9,132,372 | 9,463,942 | 1.0363 | 5,688 | 0.9978 | 0.9844 | 1.0055 | 0.0211 |
| align vs quant (reference) | pooled | 0.9686 | 0.9293 | 52,059,153 | 53,986,917 | 1.0370 | 34,094 | 0.9974 | 0.9838 | 1.0056 | 0.0219 |

Summary of ranges over the 6 samples: Spearman authors vs quant 0.9445 to 0.9491, authors vs align
0.9746 to 0.9772, align vs quant 0.9671 to 0.9701 on these 6,008 genes. On all 6,477 genes the
job reproduced the earlier align vs quant values (0.962319 to 0.966670, `validation/amb/README.md`).

Equality is not expected. The authors counted reads on CDS features with HTSeq (union mode; reads
that overlap more than one CDS gene are skipped; one count per read or pair as HTSeq counts it). We
count exons per gene: featureCounts on exons with `-p --countReadPairs` (align), or kallisto
transcript counts summed per gene (quant). I did not measure how often the RefSeq exon features
extend beyond the CDS. ncRNA genes are in our set (410 `YNC` ids, section 1). Reference annotations also differ (RefSeq GFF3 versus SGD GFF of
R64-1-1). I did not test which of these causes the 12 to 17 % higher totals on our side.

## 3. DE agreement

Authors' columns used: `AMB_only_logFC`, `AMB_only_PValue`, `AMB_only_FDR` (also present:
`AMB_only_LR`, not used). Ours: `log2FoldChange`, `pvalue`, `padj` of `deseq2.results.csv`
(DESeq2 Wald, `~ condition`, contrast `condition,AmpB,Control`, so positive = higher in AmpB).

Sign convention (`results/sign_convention.tsv`, `results/sign_convention_top10.tsv`). Two checks:

1. From the authors' own counts (6 AMB-only columns): log2((mean CPM AmpB + 0.5) / (mean CPM
   Control + 0.5)) against `AMB_only_logFC`: Spearman 0.9943 over 5,456 genes; of the 1,828 genes
   with `AMB_only_FDR` < 0.05, all 1,828 have the same sign.
2. Against our quant DESeq2 (AmpB over Control): Spearman 0.9884 over 5,405 genes; of the 1,816
   authors' FDR < 0.05 genes with a log2FC on our side, 1,816 have the same sign. The 10 genes with
   the smallest authors' FDR (YGL117W, YER069W, YDL170W, YBR147W, YER175C, YOL058W, YOR302W,
   YJL088W, YGL059W, YBR145W) all have positive values on both sides (authors 1.67 to 3.25; ours
   1.69 to 3.28).

So `AMB_only_logFC` is AMB-only over control (positive = higher with AMB), the same direction as our
contrast. No sign flip was applied.

Significance: FDR or padj < 0.05. A gene with NA padj counts as not significant (the authors'
AMB_only columns have no NA). "common" = genes in both tables. log2FC correlations use the common
genes with a non-NA log2FC on both sides. "sig, common" counts only common genes; Jaccard =
both / (x or y) on common genes. "either, same sign" = genes significant in x or y (common, log2FC
on both sides) whose two log2FC signs agree. `results/de_agreement.tsv`:

| x vs y | genes x / y | common (with log2FC both) | Spearman log2FC | Pearson log2FC | NA padj x / y | sig x / y (whole table) | sig x / y (common) | both | x only | y only | Jaccard | both, same sign | either, same sign |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| authors published vs DESeq2 quant | 5,456 / 6,477 | 5,407 (5,405) | 0.9884 | 0.9566 | 0 / 335 | 1,828 / 2,042 | 1,816 / 1,955 | 1,657 | 159 | 298 | 0.7838 | 1,657 of 1,657 | 2,114 of 2,114 |
| authors published vs DESeq2 align | 5,456 / 6,477 | 5,407 (5,398) | 0.9897 | 0.9933 | 0 / 427 | 1,828 / 2,019 | 1,816 / 1,954 | 1,672 | 144 | 282 | 0.7969 | 1,671 of 1,672 | 2,096 of 2,097 |
| DESeq2 align vs DESeq2 quant (reference), all genes | 6,477 / 6,477 | 6,477 (6,287) | 0.9265 | 0.8464 | 427 / 335 | 2,019 / 2,042 | 2,019 / 2,042 | 1,956 | 63 | 86 | 0.9292 | 1,954 of 1,956 | 2,091 of 2,101 |
| DESeq2 align vs DESeq2 quant, authors' 5,456 DE genes only | 5,407 / 5,407 | 5,407 (5,398) | 0.9964 | 0.9964 | 10 / 4 | 1,954 / 1,955 | 1,954 / 1,955 | 1,907 | 47 | 48 | 0.9525 | 1,907 of 1,907 | 1,999 of 2,000 |

The gene sets differ between rows: the authors' DE table has 5,456 genes (a subset of their 6,692
counted genes; it holds no gene with fewer than 273 reads summed over the 6 AMB samples, but it
omits genes with up to 7,217 reads, and 305 of the 1,236 omitted genes are all-zero in these 6
samples). Their filter is not stated, and it may use all 18 samples. Rows on the full 6,477 genes
include many low-count genes; the last row restricts our two results to the authors' genes for a
like-for-like reference.

## 4. Decomposition: our DESeq2 step on the authors' counts

`qrnamap de --counts $DATA/geo/prep/authors_counts_amb.tsv --samples $DATA/samples.tsv --design
'~ condition' --contrast condition,AmpB,Control --out $DATA/geo/de_authors/deseq2 --tmpdir
"$SCRATCH/geo_compare"` (job 32322903). The matrix holds all 6,692 authors' genes, header `gene` +
our sample names, plain TSV. Result: 6,692 genes, 682 NA padj, 1,959 with padj < 0.05.

(i) Method difference only (same counts: the authors'):

| x vs y | common (log2FC both) | Spearman | Pearson | sig x / y (common) | both | x only | y only | Jaccard | both, same sign | either, same sign |
|---|---|---|---|---|---|---|---|---|---|---|
| authors published vs DESeq2 on authors' counts | 5,456 (5,456) | 0.9943 | 0.9944 | 1,828 / 1,930 | 1,686 | 142 | 244 | 0.8137 | 1,686 of 1,686 | 2,072 of 2,072 |

(ii) Counting difference only (same method: our DESeq2 step):

| x vs y | gene set | common (log2FC both) | Spearman | Pearson | sig x / y (common) | both | x only | y only | Jaccard | both, same sign | either, same sign |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DESeq2 authors' counts vs DESeq2 quant | all common | 6,008 (5,884) | 0.9623 | 0.8897 | 1,941 / 2,013 | 1,874 | 67 | 139 | 0.9010 | 1,874 of 1,874 | 2,059 of 2,070 |
| DESeq2 authors' counts vs DESeq2 align | all common | 6,008 (5,875) | 0.9726 | 0.9487 | 1,941 / 1,993 | 1,891 | 50 | 102 | 0.9256 | 1,890 of 1,891 | 2,035 of 2,037 |
| DESeq2 authors' counts vs DESeq2 quant | authors' DE genes | 5,407 (5,405) | 0.9934 | 0.9580 | 1,919 / 1,955 | 1,859 | 60 | 96 | 0.9226 | 1,859 of 1,859 | 2,015 of 2,015 |
| DESeq2 authors' counts vs DESeq2 align | authors' DE genes | 5,407 (5,398) | 0.9948 | 0.9957 | 1,919 / 1,954 | 1,872 | 47 | 82 | 0.9355 | 1,871 of 1,872 | 1,998 of 1,999 |

On the authors' 5,456 DE genes, the significant-gene Jaccard is 0.8137 for (i) (method only), 0.9226
(quant) and 0.9355 (align) for (ii) (counting only), and 0.7838 (quant) and 0.7969 (align) for the
full difference (item 3). I did not test why the Pearson of log2FC is lower for quant than for
align in every comparison with quant (0.9566 vs 0.9933; 0.9580 vs 0.9957); I did not list the genes
that drive it.

## 5. Limits and caveats

- Different counting features (CDS versus exons), different counting rules (HTSeq union, ambiguous
  CDS reads skipped; featureCounts pairs; kallisto EM), different annotations (SGD GFF R64-1-1
  versus RefSeq GCF_000146045.2 GFF3) and different gene universes (6,008 common ids).
- Different DE methods. The authors' method is not stated in the GEO records. Their DE fit may use
  all 18 samples (shared dispersion or a gene filter over all samples); I did not check this.
  Our DESeq2 uses only the 6 AMB samples.
- The authors' tables: in the DE file, `LF_only_logFC` and `AMB_only_logFC` have Pearson 0.9976
  (5,363 rows without NA), median |difference| 0.0155, max 1.213. In the CPM file, max |AMB_only -
  LF_only| is 0.008 and max |control_for_AMB_only - control_for_LF_only| is 0.007; for comparison
  max |AMB_only - AMB_LF| is 2.12 and the Pearson of `AMB_only_logFC` with `AMB_and_LF_logFC` is
  0.4211 (`results/authors_table_checks.tsv`). I did not investigate this. It does not affect the
  AMB_only columns used here: their logFC matches the AMB-only columns of the authors' counts
  (Spearman 0.9943, section 3).
- Replicate-level normalization: the ratio statistics in section 2 use raw counts; library size
  differences are in them.
- Everything parsed. The 93 `NA` cells in the LF_only columns read as NA (not used).

## Defects found

No defect that changes a result. One documentation item in the CLI:

- `qrnamap de --help` prints no description for `--samples` (the line is `--samples <SAMPLES>`
  with empty help text). Reproduction: `qrnamap de --help` with the binary of `d64f917`. Other
  options of `de` have a description. Not changed (no CLI change in this task).

## Not measured

- The authors' CPM table was not compared with our values (only the column checks above).
- The AMB + LF and LF-only treatments and their controls.
- Which of the counting differences (features, UTRs, ncRNA genes, multi-mapping, ambiguity rule,
  annotation) explains the higher totals and the lower quant Pearson.
- Whether the authors' gene filter or DE method used all 18 samples.
- Whether unmatched pairs such as YDR134C (ours) and YDR133C (theirs) are the same locus under
  different names; no id translation was attempted (exact string match only).
- Shrunken log2FC (`--shrink`) and any threshold other than 0.05.
