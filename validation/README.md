# Validation on public yeast RNA-seq (Task 14)

Run on 2026-11-19 (US Pacific time) on the UCR HPCC cluster. The `qrnamap` binary was built
from commit `81175d5` (`cargo build --release -p qrnamap-cli`, pixi env on PATH).
All measured numbers in this file come from that `81175d5` build. A later commit, `f2bfb54`, changed
only the SLURM script generation of `qrnamap run --slurm` (it also added the `--tool-path` flag for
it); it fixes finding 1 below. The `f2bfb54` build was not re-run on the cluster.
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

- `prefetch`: 13.97 to 21.68 s per run. I timed it with `/usr/bin/time` in the shell. The times are
  kept in `$DATA/logs/prefetch_RUN.log` (last line). Those logs also say that prefetch fetched "SRA
  Normalized Format files with full base quality scores". `.sra` sizes 333,038,452 to 393,509,537 bytes (2.14 GB in total).
- `parallel-fastq-dump` jobs: 1:03 to 1:13 elapsed, MaxRSS 1.72 to 1.82 GB.
- Read count check: both mates of all 6 runs have exactly the ENA `read_count` (12 of 12 files).
  The `.sra` files were deleted after the check. FASTQ total: 3,428,536,290 bytes (12 `.fastq.gz`).
- Reads are up to 101 bases. The quality strings use binned NovaSeq values (`F`, `:`, `,`, `#`).
- `qrnamap fetch`: 1.75 s. I timed it with `/usr/bin/time` in the shell. No log file keeps this
  number (it was printed to the terminal only). Both files have exactly the catalog byte length:
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

Formula (`crates/qrnamap-core/src/strand.rs`): antisense fraction = antisense / (sense + antisense).
Ambiguous fragments are not in the denominator. Example, WT_1: 89,500 / (1,130 + 89,500) = 0.9875.
With the ambiguous fragments in the denominator, the value is 89,500 / (1,130 + 89,500 + 936) = 0.9774.

All 6 suggestions are `reverse`. ENA names the kit (TruSeq stranded mRNA) but not the orientation.
`reverse` is the orientation I expect for that kit (see Data). Both mode runs used `--strandedness reverse`.
strand-check job: 50 s elapsed, MaxRSS 247,472 KB.

### Per-sample results (`--samples-per-task 1`, `--threads 1`)

SLURM job ids of the reported run (the resubmission):

| mode | prepare | sample array | final |
|---|---|---|---|
| quant | 32242362 | 32242363 (tasks 0 to 5) | 32242364 |
| align | 32242365 | 32242366 (tasks 0 to 5) | 32242367 |

The first attempt failed (see finding 1). Its prepare jobs were 32242271 (quant) and 32242274
(align). They failed after 30 s because a tool was not on PATH. Their array and final jobs
(32242272, 32242273, 32242275, 32242276) never started; I cancelled them. The failed prepare jobs
did not affect the results. The second prepare job removes old outputs and writes the reference
files and the index again. Modification times (checked with `find ... -printf '%TT'`, which reads
the same mtime as `stat`), 2026-11-19:

- Result and reference files: all have an mtime of 21:47:08 or later. This is after the failed
  jobs ended (21:41:27). The earliest is `run_align/splice_sites.txt` at 21:47:08.402.
  - `run_quant/`: `transcripts.fa.zst`, `tx2gene.tsv.zst` (21:47:11), `transcripts.kidx`,
    `.prepared` (21:47:54), `SAMPLE/abundance.tsv`, `SAMPLE/run_info.json`,
    `SAMPLE.genecounts.tsv.zst` (21:50:51 to 21:52:48), `counts.tsv.zst` (21:53:17),
    `deseq2.results.csv`, `deseq2.normalized.csv`, `deseq2.vst.csv`, `deseq2.versions.txt`
    (21:53:58 to 21:53:59).
  - `run_align/`: `splice_sites.txt`, `genes.saf` (21:47:08), `index/genome.1.ht2` to
    `index/genome.8.ht2` (21:47:09 to 21:47:26), `.prepared` (21:47:26),
    `SAMPLE/featurecounts.summary`, `SAMPLE.genecounts.tsv.zst` (21:55:03 to 21:58:44),
    `counts.tsv.zst` (21:58:56), `deseq2.*` (21:59:22 to 21:59:23).
- Exceptions (earlier than 21:47):
  - The 4 generated scripts in each `slurm/` directory (`prepare.sbatch`, `samples.sbatch`,
    `final.sbatch`, `submit.sh`): 21:40:31. They were written once by `qrnamap run --slurm`
    and not rewritten for the resubmission. Both attempts used the same scripts.
  - The first-attempt prepare logs `qrnamap-prepare_32242271.out` and
    `qrnamap-prepare_32242274.out`: 21:41:27.
  - The 0-byte logs of the reported prepare jobs, `qrnamap-prepare_32242362.out` and
    `qrnamap-prepare_32242365.out`: 21:45:19.

What the log directories contain (`$DATA/run_MODE/slurm/logs/`):

- `qrnamap-prepare_32242271.out` (quant) and `qrnamap-prepare_32242274.out` (align): the error line
  of the failed first attempt. These are the only logs from the first attempt.
- `qrnamap-prepare_32242362.out` and `qrnamap-prepare_32242365.out`: empty (0 bytes). The prepare
  step of the reported run printed nothing.
- `qrnamap-samples_JOB_N.out` (6 per mode) and `qrnamap-final_JOB.out`: the reported run.

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
- Align assigned fraction = featureCounts `Assigned` / read pairs processed by HISAT2 (the
  `reads_processed` value in the CLI log, from the HISAT2 line `N reads; of these:`). The
  denominator is not the sum of the featureCounts summary. That sum is larger, because featureCounts
  counts each alignment of a multi-mapping pair. Example, WT_1: Assigned 5,031,627; pairs processed
  5,618,967; sum of `featurecounts.summary` 6,279,913. Assigned / pairs = 0.8955 (reported).
  Assigned / summary sum = 0.8012. Assigned / summary sum for all 6 samples: 0.7923 to 0.8012.
- Quant task 0 (WT_1) ran on node c07 with 1 CPU allocated. Tasks 1 to 5 ran on r42 with 2 CPUs
  allocated (SLURM gave 2; the scripts ask for 1, and kallisto ran with `-t 1`). WT_1 took 264 s;
  the others took 147 to 173 s. I did not test whether the node or the sample caused the difference.
  All align tasks ran on c07 at the same time, each with 1 CPU.
- featureCounts summary for WT_1 (align): Assigned 5,031,627; Unassigned_MultiMapping 1,010,770;
  Unassigned_NoFeatures 168,715; Unassigned_Ambiguity 44,381; Unassigned_Unmapped 24,420.
- Other jobs: quant prepare (transcripts and kallisto index) 2:36; align prepare (SAF, splice
  sites, HISAT2 index) 2:08; quant final (matrix and DESeq2) 0:44; align final 0:29.
- No `warning:` line appeared in any job log.

### Bundling samples per task (derived from the times above)

The run used `--samples-per-task 1`. The values 20 and 8 below are a recommendation derived from
the measured times. They are not a tested setting. A run with the recommended values was not
done; the run used 1 sample per task.

Measured per-sample time, quant: mean 174.3 s, max 264 s. The quant mean includes both node types
(1 task on c07, 5 tasks on r42). The max (264 s) is the WT_1 task on c07, so the max-based bound is
also shown. Align: mean 583.3 s, max 652 s.
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
   32242271 and 32242274. I did not change the CLI during this run. I cancelled the dependent
   jobs and resubmitted the same scripts with the env on PATH (stage `submit`), which is the
   same requirement as in local mode.

   **Status: fixed in f2bfb54 after this run.** Commit f2bfb54 ("fix(cli): SLURM scripts bake
   tool paths, check tools, kill on invalid dependency (R28)") makes `run --slurm` write the tool
   directories into `export PATH=` in every `.sbatch` file. `prepare.sbatch` checks the tools and
   exits with a list of the missing ones. `submit.sh` passes `--kill-on-invalid-dep=yes` on the
   dependent jobs. The fix was verified by the project's unit tests with a stub `sbatch` and by
   `sbatch --test-only` (option parsing only). I did not run those checks myself in this task.
   The fix was NOT re-run on the cluster. The claim that the fixed scripts work on real jobs
   without the submitter's PATH is untested here.
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
