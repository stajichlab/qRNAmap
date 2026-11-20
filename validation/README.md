# Validation on public yeast RNA-seq (Task 14)

Run on 2026-11-19 (US Pacific time) on the UCR HPCC cluster. The `qrnamap` binary was built
from commit `81175d5` (`cargo build --release -p qrnamap-cli`, pixi env on PATH).
Data lives outside the repository in `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/`.
This directory holds the scripts and small summary tables only.

Files:

- `run_yeast.sh`: the driver. Each stage is one command (see the header of the script).
- `compare_gffread.py`: compares two transcript FASTA files by id and sequence (reads plain, `.gz`, `.zst`).
- `compare_modes.R`: count and DE agreement between the two modes.
- `results/`: the tables quoted below (`ena_runs.tsv`, `parity_summary.txt`, `strand_check.tsv`,
  `per_sample.tsv`, `mode_agreement.tsv`, `de_overlap.tsv`).

## Data

Study: ENA `PRJNA674653` (SRA `SRP291191`, GEO `GSE160861`), "Disturbed intramitochondrial
phosphatidic acid transport impairs cellular stress signaling" (PubMed 33497623). Strain S288c,
Illumina NovaSeq 6000, paired-end, `library_strategy` RNA-Seq, `library_source` TRANSCRIPTOMIC,
`library_selection` cDNA. The study has 12 runs: WT and ups1 deletion, each untreated and
tunicamycin-treated, 3 biological replicates each. I used the 6 untreated runs (WT versus ups1 KO).

How I chose it: ENA portal API `search` with `tax_eq(4932)`, RNA-Seq, TRANSCRIPTOMIC, PAIRED,
ILLUMINA, `read_count` between 3,000,000 and 8,000,000 (3,915 runs). Then I looked at studies with
6 to 12 runs in that range and an S288C-derived strain. I did not inspect every study. Of the
studies I inspected, PRJNA674653 is the only one with a clean 2 x 3 design in which all 6 runs are
under 8 million pairs and the strain field is `S288c`. Other candidates had a non-S288C strain (W303-type anchor-away strain HHY221 in PRJNA880888,
CEN.PK in PRJNA491917), or a replicate above 8 million pairs (BY4741, PRJNA694910).

Condition labels come from the ENA `sample_title` (`WT-___-____` = WT, untreated; `KO-___-ups1` =
ups1 KO, untreated). ENA values (`filereport`, `fields=read_count,base_count,library_layout`):

| run | sample_title | condition | read_count (pairs) | base_count | layout |
|---|---|---|---|---|---|
| SRR12982935 | S_001-F_AKIN-L____WT-___-____-REP_1 | WT_1 | 5,618,967 | 1,111,744,150 | PAIRED |
| SRR12982936 | S_002-F_AKIN-L____WT-___-____-REP_2 | WT_2 | 5,337,611 | 1,058,146,138 | PAIRED |
| SRR12982937 | S_003-F_AKIN-L____WT-___-____-REP_3 | WT_3 | 6,088,168 | 1,209,210,943 | PAIRED |
| SRR12982938 | S_004-F_AKIN-L____KO-___-ups1-REP_1 | ups1KO_1 | 5,156,073 | 1,025,112,924 | PAIRED |
| SRR12982939 | S_005-F_AKIN-L____KO-___-ups1-REP_2 | ups1KO_2 | 5,444,191 | 1,082,887,079 | PAIRED |
| SRR12982940 | S_006-F_AKIN-L____KO-___-ups1-REP_3 | ups1KO_3 | 5,311,644 | 1,055,850,537 | PAIRED |

Library kit (ENA `library_construction_protocol`): "Illumina TruSeq mRNA stranded sample
preparation Kit", poly-A selection, random-primed first strand, second strand with DNA Polymerase I
and RNase H. ENA does not state the read orientation. The orientation expected for this kit
(read 1 antisense, `reverse`) is my knowledge of the kit, not a value from ENA.

Reference: NCBI RefSeq `GCF_000146045.2` (R64), annotation SGD R64-5-1 (catalog entry).

## Commands

All stages: `bash /bigdata/stajichlab/jstajich/projects/qRNAmap/.claude/worktrees/cli-core/validation/run_yeast.sh STAGE`.
Compute stages ran as `sbatch` jobs on the `short` partition. `epyc` had a long queue.

1. Reads: `prefetch RUN -O $DATA/sra --max-size 30G` on the session node (sratoolkit 3.2.0 module),
   then one SLURM job per run (8 CPUs, 16 GB): `parallel-fastq-dump --sra-id $DATA/sra/RUN/RUN.sra
   --threads 8 --outdir $DATA/fastq --split-files --gzip --tmpdir "${SCRATCH:?}"`
   (parallel-fastq-dump 0.6.7 module). The job counts the reads of both mates and fails on a
   mismatch with the ENA `read_count`. It deletes the `.sra` after the check passes.
   Stages `prefetch`, `submit-dump` (`dump RUN PAIRS` is the job body).
2. Genome: `qrnamap fetch --source ncbi --id GCF_000146045.2 --out $DATA/genome` (stage `genome`).
3. Parity (stage `parity`, SLURM): `qrnamap transcripts --reference GENOME.fna.gz --gff GENOME.gff.gz
   --out qrnamap_transcripts.fa.zst --tx2gene tx2gene.tsv.zst`;
   `gffread -w gffread_transcripts.fa -g genome.fa genes.gff3` (plain copies of the same two files
   in `$SCRATCH`); `compare_gffread.py --ours ... --gffread ... --gff ... --out parity_diff.tsv`.
4. Sample sheet (stage `sheet`): columns `sample fastq1 fastq2 condition replicate`.
5. Strandedness (stage `strand`, SLURM, 4 CPUs): `qrnamap strand-check --samples samples.tsv
   --reference ... --gff ... --out $DATA/strand --threads 4` (default `--pairs 100000`).
6. Both modes (stages `slurm MODE 1 reverse` then `submit MODE`):
   `qrnamap run --samples samples.tsv --reference ... --gff ... --mode MODE --out $DATA/run_MODE
   --strandedness reverse --threads 1 --contrast condition,ups1KO,WT --slurm --samples-per-task 1
   --partition short --mem-gb 16 --time 02:00:00`, then `bash $DATA/run_MODE/slurm/submit.sh`
   with the pixi env on PATH.
7. Agreement (stage `compare`, SLURM): `Rscript compare_modes.R $DATA/run_quant $DATA/run_align results`.

Tool versions (pixi env `default`): kallisto 0.52.0, HISAT2 2.2.2 (`hisat2-align-s`),
featureCounts 2.1.1, gffread 0.12.9, R 4.5.3, DESeq2 1.50.2, zstd 1.5.7, Python 3.15.2.
Cluster modules: sratoolkit 3.2.0, parallel-fastq-dump 0.6.7.

## Measured

All numbers below come from the runs of 2026-11-19 described above.

### Download

- `prefetch`: 13.97 to 21.68 s per run. `.sra` sizes 333,038,452 to 393,509,537 bytes (2.14 GB in total).
- `parallel-fastq-dump` jobs: 1:03 to 1:13 elapsed, MaxRSS 1.72 to 1.82 GB.
- Read count check: both mates of all 6 runs have exactly the ENA `read_count` (12 of 12 files).
  The `.sra` files were deleted after the check. FASTQ total: 3,428,536,290 bytes (12 `.fastq.gz`).
- Reads are up to 101 bases. The quality strings use binned NovaSeq values (`F`, `:`, `,`, `#`).
- `qrnamap fetch`: 1.75 s. Both files have exactly the catalog byte length:
  `GCF_000146045.2_R64_genomic.fna.gz` 3,843,460 bytes, `GCF_000146045.2_R64_genomic.gff.gz`
  2,203,839 bytes.

### Transcript parity with gffread

`results/parity_summary.txt`: our FASTA has 6,478 records, gffread has 6,478 records,
6,478 are identical (by id and upper-cased sequence). Missing on our side: 0. Missing on the
gffread side: 0. Different sequence: 0. Duplicate ids: 0 on each side. There is no difference to
explain. The tx2gene table has 6,478 transcripts and 6,477 genes (gene `Q0158` has 2 transcripts).
Parity job: 8 s elapsed.

### Strandedness

`results/strand_check.tsv` (first 100,000 pairs per sample):

| sample | sense | antisense | ambiguous | antisense fraction | suggested |
|---|---|---|---|---|---|
| WT_1 | 1,130 | 89,500 | 936 | 0.9875 | reverse |
| WT_2 | 1,086 | 89,636 | 737 | 0.9880 | reverse |
| WT_3 | 1,261 | 89,803 | 668 | 0.9862 | reverse |
| ups1KO_1 | 1,305 | 88,724 | 1,217 | 0.9855 | reverse |
| ups1KO_2 | 1,194 | 88,758 | 1,475 | 0.9867 | reverse |
| ups1KO_3 | 1,194 | 88,818 | 1,464 | 0.9867 | reverse |

All 6 suggestions are `reverse`. ENA names the kit (TruSeq stranded mRNA) but not the orientation.
`reverse` is the orientation I expect for that kit (see Data). Both mode runs used `--strandedness reverse`.
strand-check job: 50 s elapsed, MaxRSS 247,472 KB.

### Per-sample results (`--samples-per-task 1`, `--threads 1`)

`results/per_sample.tsv`. Elapsed and MaxRSS come from `sacct` (the `.batch` step of each array
task). One array task holds one sample. The time includes reading the FASTQ, the tool runs and
writing the outputs.

| sample | quant: pseudoaligned % | quant: elapsed | align: HISAT2 overall alignment % | align: assigned fraction | align: elapsed | align: MaxRSS (KB) |
|---|---|---|---|---|---|---|
| WT_1 | 93.60 | 4:24 | 96.84 | 0.8955 | 10:16 | 7,295,044 |
| WT_2 | 94.30 | 2:27 | 97.23 | 0.8966 | 10:02 | 6,968,312 |
| WT_3 | 94.60 | 2:53 | 97.53 | 0.8979 | 10:52 | 7,953,212 |
| ups1KO_1 | 93.10 | 2:32 | 96.81 | 0.8878 | 9:52 | 6,727,536 |
| ups1KO_2 | 93.20 | 2:35 | 97.27 | 0.8875 | 10:07 | 7,105,324 |
| ups1KO_3 | 93.20 | 2:35 | 96.65 | 0.8879 | 7:11 | 6,917,516 |

- In quant mode the assigned fraction is the pseudoaligned fraction (kallisto reports
  `p_pseudoaligned` with one decimal). Quant MaxRSS: 107,420 to 111,888 KB.
- Quant task 0 ran on node c07; tasks 1 to 5 ran on r42. All align tasks ran on c07 at the same time.
- featureCounts summary for WT_1 (align): Assigned 5,031,627; Unassigned_MultiMapping 1,010,770;
  Unassigned_NoFeatures 168,715; Unassigned_Ambiguity 44,381; Unassigned_Unmapped 24,420.
- Other jobs: quant prepare (transcripts and kallisto index) 2:36; align prepare (SAF, splice
  sites, HISAT2 index) 2:08; quant final (matrix and DESeq2) 0:44; align final 0:29.
- No `warning:` line appeared in any job log.

### Bundling samples per task (derived from the times above)

Measured per-sample time, quant: mean 174.3 s, max 264 s. Align: mean 583.3 s, max 652 s.
For a task of 1 to 1.5 hours:

- quant: `--samples-per-task 20` (20 x 174.3 s = 58 min at the mean; 20 x 264 s = 88 min at the max).
- align: `--samples-per-task 8` (8 x 583.3 s = 78 min at the mean; 8 x 652 s = 87 min at the max).

These values are for about 5.2 to 6.1 million pairs per sample, `--threads 1`, on nodes c07 and r42.
With 6 samples, either value gives one task for this data set. I did not run a bundled run.
Time per million pairs: quant 27.5 to 47.0 s, align 81.1 to 114.8 s.

### Agreement between modes

Both gene matrices have 6,477 rows. The row names are identical and in the same order, and the
column names are identical (checked in `compare_modes.R` with `identical()`).

Spearman correlation of gene counts between quant and align, per sample, across all 6,477 genes.
Computed with R `cor(method = "spearman")` (ties get average ranks). `results/mode_agreement.tsv`:

| sample | Spearman, all genes | Spearman, genes nonzero in either mode | Pearson of log2(count+1) | total quant | total align |
|---|---|---|---|---|---|
| WT_1 | 0.972407 | 0.968034 (6,086 genes) | 0.964266 | 5,257,657 | 5,031,627 |
| WT_2 | 0.973373 | 0.969306 (6,096) | 0.965242 | 5,032,232 | 4,785,858 |
| WT_3 | 0.973390 | 0.969712 (6,138) | 0.965536 | 5,759,433 | 5,466,754 |
| ups1KO_1 | 0.971497 | 0.965747 (6,013) | 0.963468 | 4,802,083 | 4,577,580 |
| ups1KO_2 | 0.971456 | 0.966122 (6,049) | 0.964074 | 5,076,600 | 4,831,553 |
| ups1KO_3 | 0.972502 | 0.966902 (6,018) | 0.964710 | 4,951,175 | 4,716,187 |
| all 38,862 values pooled | 0.972466 | 0.967705 (36,400) | 0.964567 | 30,879,180 | 29,409,559 |

Per gene across the 6 samples (genes with nonzero variance in both modes, 6,114 genes): median
Spearman 1. With 6 values per gene this number has little resolution.

DESeq2 (`~ condition`, Wald, contrast `condition,ups1KO,WT`, no shrinkage), genes with
`padj < 0.05` (`results/de_overlap.tsv`):

| measure | value |
|---|---|
| genes in each result | 6,477 (quant), 6,477 (align) |
| padj NA | 496 (quant), 400 (align) |
| padj < 0.05, quant | 2,648 |
| padj < 0.05, align | 2,632 |
| in both | 2,571 |
| quant only | 77 |
| align only | 61 |
| Jaccard (both / either) | 0.949059 |
| in both with the same log2FC sign | 2,571 of 2,571 |
| Spearman of log2FC, genes with a log2FC in both | 0.955239 |

The plan has no pass or fail criterion for these numbers. I do not call them good or bad.

## Findings and defects found

1. **The generated SLURM scripts depend on the PATH of the shell that runs `submit.sh`.**
   They call kallisto, hisat2-build, featureCounts and Rscript by name. They do not set or check
   PATH, and neither the scripts nor the `run --slurm` message says that the tools must be on
   PATH at submit time. Reproduction: in a shell without the pixi env on PATH, run
   `qrnamap run ... --mode quant --slurm` (the qrnamap process itself had the env on PATH) and then
   `bash OUT/slurm/submit.sh`. The prepare job fails after 30 s with
   `error: cannot start 'kallisto': No such file or directory (os error 2) (is it on PATH?)`
   (align mode: `cannot start 'hisat2-build'`). The array and final jobs then stay PENDING with
   reason `DependencyNeverSatisfied` until they are cancelled by hand. Observed with jobs
   32242271 and 32242274. I did not change the CLI. I resubmitted with the env on PATH
   (stage `submit`), which is the same requirement as in local mode. Possible fixes (not done):
   write the PATH of `qrnamap run --slurm` into the scripts, or check the tools when the scripts
   are written.
2. **Align tasks have a MaxRSS of 6.7 to 8.0 GB with `hisat2 -p 1`** (quant tasks: about 0.1 GB).
   The SAM goes to a file in `$SCRATCH`, not to memory. I did not find the cause. The cluster's
   memory accounting may count the page cache of the files written in the job; this is UNVERIFIED.
   The 16 GB request was enough.

## Not measured, or not done

- No ground truth exists for this data. The numbers show agreement between the two modes and
  with gffread. They do not show accuracy.
- Local mode (`qrnamap run` without `--slurm`) was not run on this data.
- No bundled run (`--samples-per-task` above 1). The values above are a recommendation only.
- Memory with bundled samples, and run time at other read depths or thread counts.
- The cause of the align MaxRSS (finding 2).
- strand-check used only the first 100,000 pairs per sample (the default).
- `--threads` above 1 and `--keep-sam` were not used.
- DE with LFC shrinkage, and the tunicamycin-treated runs of the study, were not used.
- No comparison against the published results of the study.
- The read orientation for the TruSeq stranded mRNA kit is not stated in ENA. The match with
  `reverse` relies on my knowledge of the kit.
