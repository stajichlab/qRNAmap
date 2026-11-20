# AMB experiments on the cluster (E0 to E3)

Run on 2026-11-20 (US Pacific time) on the UCR HPCC cluster (SLURM 24.11.1). The `qrnamap` binary
was built from commit `ff611ae` (branch `harden-cli`, created from `main` at `ff611ae`) with
`cargo build --release -p qrnamap-cli`, the pixi env on PATH, and
`CARGO_TARGET_DIR=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-harden`. No CLI code
was changed in this task. All numbers in this file come from the runs described here.

Data and run outputs are outside the repository, in
`/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb/` (called `$DATA` below). This directory
holds the scripts and small tables only.

Files:

- `e0_integrity.sh`: E0 job body (`gzip -t` and record count of each FASTQ file).
- `run_amb.sh`: the driver of E1 to E3. Each stage is one command (see the header of the script).
- `collect_amb.py`: the per-sample table of E2 (`results/per_sample.tsv`).
- `compare_modes.R`: a copy of `validation/compare_modes.R` (same logic): count and DE agreement.
- `results/`: `ena_runs.tsv`, `e0_fastq_counts.tsv`, `strand_check.tsv`, `per_sample.tsv`,
  `mode_agreement.tsv`, `de_overlap.tsv`, `e3_memory.tsv`, `jobs.tsv` (sacct of every job).

## Data

GEO `GSE80357`, BioProject `PRJNA318684`, *Saccharomyces cerevisiae* S288C, RPMI-1640, 30 C.
AMB = amphotericin B 0.03 ug/mL for 60 min. The matched controls get no drug. Illumina HiSeq 2000,
paired-end. FASTQ files (read-only, not changed by this task):
`/bigdata/gen220/shared/data-examples/rnaseq/kallisto/data/<RUN>_1.fastq.gz` and `_2.fastq.gz`.
Sample sheet source: `samples_AMB.tsv` in that directory. The top-level `samples.tsv` there was
NOT used (its controls SRR3396416 to SRR3396418 are "Matched control for LF only", single-end
NextSeq).

ENA values (`filereport`, `results/ena_runs.tsv`):

| sample | run | ENA sample_title | read_count (pairs) | base_count |
|---|---|---|---|---|
| Control_1 | SRR3396385 | Matched control for AMB only rep 1 | 9,625,471 | 1,944,345,142 |
| Control_2 | SRR3396386 | Matched control for AMB only rep 2 | 10,302,142 | 2,081,032,684 |
| Control_3 | SRR3396387 | Matched control for AMB only rep 3 | 9,997,148 | 2,019,423,896 |
| AmpB_1 | SRR3396381 | AMB only rep 1 | 8,405,534 | 1,697,917,868 |
| AmpB_2 | SRR3396382 | AMB only rep 2 | 9,755,935 | 1,970,698,870 |
| AmpB_3 | SRR3396384 | AMB only rep 3 | 10,175,785 | 2,055,508,570 |

All six: `PAIRED`, `Illumina HiSeq 2000`, `RNA-Seq`, `cDNA`, `TRANSCRIPTOMIC`.

Library kit. The ENA `library_construction_protocol` (same text on the runs) says: "For the AMB
only treatment, the AMB+LF treatment and their corresponding matching controls, the TruSeq RNA
Library prep kit v2 (Illumina) with a poly-A pulldown was used". For the LF-only runs it names a
different kit, "TruSeq Stranded mRNA-seq prep kit". The record does not say whether the TruSeq RNA
v2 kit is stranded. That it is not stranded is my knowledge of the kit, not a value from the record.
The authors counted with HTSeq `-s no` (unstranded).

Sample sheet `$DATA/samples.tsv`: columns `sample fastq1 fastq2 condition replicate`; `condition`
is `Control` or `AmpB`; `replicate` is `r1` to `r3` (from `Rep1` to `Rep3` of `samples_AMB.tsv`).
FASTQ paths are the absolute paths of the shared files.

Reference (reused, read-only): NCBI RefSeq `GCF_000146045.2` (R64),
`/bigdata/stajichlab/jstajich/projects/qRNAmap_data/genome/GCF_000146045.2_R64_genomic.fna.gz` and
`.gff.gz`.

Tool versions (pixi env `default` of this worktree, created with `pixi install`): kallisto 0.52.0,
HISAT2 2.2.2 (`hisat2-align-s`), featureCounts 2.1.1, R 4.5.3, DESeq2 1.50.2, zstd 1.5.7,
Python 3.15.2. SLURM 24.11.1. Kernel 4.18.0-553.184.1.el8_10.

## E0. Integrity of the input

Command: `sbatch -p short -c 2 --mem 2G --wrap "bash .../validation/amb/e0_integrity.sh SHARED_DIR $DATA/e0/fastq_counts.tsv"`
(job 32307800, node c27, 11:03 elapsed). Records = `zcat | wc -l` / 4. ENA `read_count` from
`https://www.ebi.ac.uk/ena/portal/api/filereport?accession=RUN&result=read_run&fields=read_count`.

| file | gzip -t exit | records | ENA read_count | match |
|---|---|---|---|---|
| SRR3396385_1 | 0 | 9,625,471 | 9,625,471 | yes |
| SRR3396385_2 | 0 | 9,625,471 | 9,625,471 | yes |
| SRR3396386_1 | 0 | 10,302,142 | 10,302,142 | yes |
| SRR3396386_2 | 0 | 10,302,142 | 10,302,142 | yes |
| SRR3396387_1 | 0 | 9,997,148 | 9,997,148 | yes |
| SRR3396387_2 | 0 | 9,997,148 | 9,997,148 | yes |
| SRR3396381_1 | 0 | 8,405,534 | 8,405,534 | yes |
| SRR3396381_2 | 0 | 8,405,534 | 8,405,534 | yes |
| SRR3396382_1 | 0 | 9,755,935 | 9,755,935 | yes |
| SRR3396382_2 | 0 | 9,755,935 | 9,755,935 | yes |
| SRR3396384_1 | 0 | 10,175,785 | 10,175,785 | yes |
| SRR3396384_2 | 0 | 10,175,785 | 10,175,785 | yes |

12 of 12 files pass `gzip -t` and have exactly the ENA `read_count`. Table: `results/e0_fastq_counts.tsv`.

## E1. D1 on real jobs, failure path

Generation (pixi env on PATH of the generating shell only): `run_amb.sh e1-gen`, that is
`qrnamap run --samples $DATA/e1/samples_missing.tsv ... --mode quant --out $DATA/e1/run_quant
--strandedness unstranded --contrast condition,AmpB,Control --slurm --samples-per-task 6
--partition short --mem-gb 4 --time 00:10:00`. The sheet points at files under
`$DATA/e1/does_not_exist/`, which does not exist. Generation did not check the files and wrote the
scripts. `prepare.sbatch` has `export PATH='<worktree>/.pixi/envs/default/bin':"$PATH"` and a
`command -v` check of `kallisto` and `Rscript`. `submit.sh` passes `--kill-on-invalid-dep=yes` to
the array and final jobs.

Two submissions of the same scripts:

| attempt | submit shell | prepare | array | final |
|---|---|---|---|---|
| E1a | `env -i HOME USER PATH=<slurm bin>:/usr/bin:/bin bash submit.sh` (`run_amb.sh submit-clean`) | 32308072 | 32308073 | 32308074 |
| E1b | `env -i HOME USER PATH=/usr/bin:/bin bash -l -c 'bash submit.sh'` (`run_amb.sh submit-login`) | 32308173 | 32308174 | 32308175 |

sacct after the failure (`sacct -j ... --format=JobID,JobName,State,ExitCode,Elapsed,Start,End,Reason`):

| job | E1a state | E1a end | E1b state | E1b end |
|---|---|---|---|---|
| prepare | FAILED 1:0, 00:00:01 | 09:45:34 | FAILED 1:0, 00:00:01 | 09:48:40 |
| samples array (task 0) | CANCELLED, never started | 09:45:46 | CANCELLED, never started | 09:48:52 |
| final | CANCELLED, never started | 09:45:46 | CANCELLED, never started | 09:48:52 |

- `scontrol show job` of the array and final jobs (both attempts): `JobState=CANCELLED
  Reason=DependencyNeverSatisfied`, `KillOInInvalidDependent=Yes`; dependency
  `afterok:<prepare>(failed)` for the array, `afterok:<array>_*(failed)` for the final job.
- SLURM cancelled both dependent jobs 12 s after the prepare job ended, in both attempts.
  `squeue -j <array>,<final>` about 1 minute after the failure (E1b, 09:50:05) listed no job.
  No job stayed pending. I cancelled nothing by hand.
- E1b prepare log: `error: sample 'Control_1': input error: $DATA/e1/does_not_exist/SRR3396385_1.fastq.gz:
  No such file or directory (os error 2)`. This is the up-front FASTQ check, as intended.
- E1a prepare log: `/var/spool/slurmd/job32308072/slurm_script: line 17: SCRATCH: SCRATCH is not
  set: this job needs node-local scratch`. The job failed at the `${SCRATCH:?}` line, before the
  FASTQ check. See defect D3. A probe job (32308141, `--wrap 'echo ${SCRATCH:-UNSET}'`, submitted
  from `env -i` with `LOADEDMODULES=workspace/scratch`) also printed `SCRATCH=UNSET`.

Result: with `--kill-on-invalid-dep=yes`, a failed prepare job made SLURM cancel the array and final
jobs (2 of 2 attempts). The flag works on this cluster.

## E2. D1 on real jobs, positive path, and the AMB data in both modes

### Strandedness

`run_amb.sh strand` as an sbatch job (job 32308227, node r46, 4 CPUs, 0:33 elapsed, MaxRSS
250,508 KB): `qrnamap strand-check --samples $DATA/samples.tsv --reference ... --gff ... --out
$DATA/strand --threads 4` (default `--pairs 100000`). `results/strand_check.tsv`:

| sample | sense | antisense | ambiguous | antisense fraction | suggested |
|---|---|---|---|---|---|
| Control_1 | 45,434 | 45,217 | 628 | 0.4988 | unstranded |
| Control_2 | 46,474 | 45,729 | 784 | 0.4960 | unstranded |
| Control_3 | 44,969 | 44,265 | 749 | 0.4961 | unstranded |
| AmpB_1 | 46,282 | 45,590 | 751 | 0.4962 | unstranded |
| AmpB_2 | 45,543 | 45,129 | 760 | 0.4977 | unstranded |
| AmpB_3 | 46,322 | 45,422 | 792 | 0.4951 | unstranded |

Antisense fraction = antisense / (sense + antisense) (ambiguous not in the denominator). Example
Control_1: 45,217 / (45,434 + 45,217) = 0.4988. Each sample also printed the warning "neither
strand holds more than 0.8 of the stranded-assigned fragments". Both mode runs used
`--strandedness unstranded`. This agrees with the authors' `-s no`. The kit evidence is in Data.

### Runs

Generation (pixi env on PATH of the generating shell only):

- quant: `run_amb.sh slurm quant 6 unstranded 8 02:00:00`
- align: `run_amb.sh slurm align 3 unstranded 32 02:00:00`

that is `qrnamap run --samples $DATA/samples.tsv --reference ... --gff ... --mode MODE --out
$DATA/run_MODE --strandedness unstranded --threads 1 --contrast condition,AmpB,Control --slurm
--samples-per-task K --partition short --mem-gb M --time 02:00:00`. The design is the default
`~ condition`. `--samples-per-task 6` (quant) gives 1 array task; `3` (align) gives 2 tasks.

`--mem-gb` and `--time` are placeholders, set from the yeast measurements with a margin: quant 8 GB
(yeast quant MaxRSS about 0.11 GB); align 32 GB (yeast 1,329 to 1,338 bytes per pair x about 10
million pairs = about 13.4 GB per sample); `--time 02:00:00` (the `short` limit; yeast maxima x
the AMB pair counts gave about 46 min for the quant task and about 57 min per align task).

Submission: `run_amb.sh submit-login $DATA/run_MODE`, a fresh login shell (`env -i HOME USER
PATH=/usr/bin:/bin bash -l`). In that shell `kallisto`, `hisat2-align-s`, `hisat2-build-s` and
`featureCounts` were NOT found, and `Rscript` was the cluster module R 4.5.2
(`/opt/linux/rocky/8.x/x86_64/pkgs/R/4.5.2/bin/Rscript`), not the pixi R. I did not use the
`env -i` shell of E1a here, because jobs submitted from it get no `$SCRATCH` (D3).

Baked PATH result: all 6 E2 jobs COMPLETED with exit 0:0. The prepare jobs ran `kallisto index` and
`hisat2-build`; the array tasks ran kallisto, `hisat2-align-s` and featureCounts; the final jobs ran
DESeq2. `deseq2.versions.txt` of both runs says `R version 4.5.3` and `DESeq2 1.50.2`: the pixi R,
not the R 4.5.2 of the submit shell. So the baked `export PATH` worked on real jobs.

| mode | prepare | sample array | final |
|---|---|---|---|
| quant | 32308256: 1:44, MaxRSS 182,448 KB, r43 | 32308257 (task 0, 6 samples): 23:29, MaxRSS 118,460 KB, r41 | 32308258: 0:29, MaxRSS 870,664 KB, x05 |
| align | 32308259: 1:29, MaxRSS 157,084 KB, r43 | 32308260 task 0 (Control_1..3): 20:51, MaxRSS 13,070,508 KB, r46; task 1 (AmpB_1..3): 22:13, MaxRSS 12,895,776 KB, r46 | 32308261: 0:35, MaxRSS 899,780 KB, r44 |

Elapsed and MaxRSS are from `sacct` (the `.batch` step). Every job asked for 1 CPU
(`--cpus-per-task=1`); SLURM allocated 2 (AllocCPUS 2). The tools ran with 1 thread. Node r46 also
ran the three E3 jobs (and both align tasks) at the same time, from 09:52 to 10:00. `results/jobs.tsv`
has every job.

Throughput from the task times: quant 1,409 s for 58,262,015 pairs = 24.2 s per million pairs.
Align task 0: 1,251 s for 29,924,761 pairs = 41.8 s per million pairs; task 1: 1,333 s for
28,337,254 pairs = 47.0 s per million pairs.

### Per sample

`results/per_sample.tsv` (`collect_amb.py`). Pseudoaligned % is kallisto `p_pseudoaligned` in
`run_info.json`. HISAT2 overall alignment % and pairs processed are from the CLI log line of the
array task. Assigned fraction (align) = featureCounts `Assigned` / read pairs processed by HISAT2.
Seconds per sample are derived from the mtimes of `SAMPLE.genecounts.tsv.zst`: the difference to
the previous sample of the same task, or to the task's sacct Start for the first sample of a task.
They are not timer values.

| sample | pairs | quant: pseudoaligned % | quant: s | align: HISAT2 overall % | align: Assigned | align: assigned fraction | align: s |
|---|---|---|---|---|---|---|---|
| Control_1 | 9,625,471 | 93.6 | 258.2 | 95.77 | 8,594,998 | 0.8929 | 341.1 |
| Control_2 | 10,302,142 | 94.6 | 254.8 | 97.06 | 9,361,327 | 0.9087 | 528.3 |
| Control_3 | 9,997,148 | 94.7 | 237.9 | 96.67 | 8,799,445 | 0.8802 | 379.7 |
| AmpB_1 | 8,405,534 | 94.0 | 191.1 | 96.95 | 7,610,562 | 0.9054 | 296.1 |
| AmpB_2 | 9,755,935 | 93.6 | 230.1 | 95.49 | 8,714,096 | 0.8932 | 508.5 |
| AmpB_3 | 10,175,785 | 94.1 | 237.1 | 96.33 | 9,173,284 | 0.9015 | 526.9 |

- kallisto `n_processed` equals the ENA pair count for all 6 samples, and so does the HISAT2 pair count.
- featureCounts summary sum (counts each alignment of a multi-mapping pair) for Control_1:
  10,328,187; Unassigned_MultiMapping 1,117,148; Unassigned_NoFeatures 351,479;
  Unassigned_Ambiguity 65,177; Unassigned_Unmapped 199,385. All samples are in the table file.
- Per-sample seconds, range: quant 191.1 to 258.2 s (22.7 to 26.8 s per million pairs); align
  296.1 to 528.3 s (35.2 to 52.1 s per million pairs). In both align tasks the first sample took
  296 to 341 s and the later ones 380 to 528 s. I did not test why.
- No line starting with `warning:` appeared in the array or final logs.

### Agreement between modes

`run_amb.sh compare` as an sbatch job (32309155, c19, 0:24): `Rscript compare_modes.R $DATA/run_quant
$DATA/run_align results`.

Both matrices have 6,477 gene rows. Row names are identical and in the same order; column names
are identical (`identical()` in `compare_modes.R`).

Spearman correlation of gene counts between quant and align, per sample over all 6,477 genes,
computed with R `cor(method = "spearman")` (ties get average ranks). Pooled = all 38,862 values
(6 samples x 6,477 genes) in one vector per mode. `results/mode_agreement.tsv`:

| sample | Spearman, all genes | Spearman, genes nonzero in either mode | Pearson of log2(count+1) | total quant | total align |
|---|---|---|---|---|---|
| Control_1 | 0.963560 | 0.962156 (6,376) | 0.928885 | 9,005,877 | 8,594,998 |
| Control_2 | 0.964354 | 0.962828 (6,367) | 0.937710 | 9,744,045 | 9,361,327 |
| Control_3 | 0.962319 | 0.960342 (6,345) | 0.934092 | 9,462,768 | 8,799,445 |
| AmpB_1 | 0.966670 | 0.965182 (6,363) | 0.940109 | 7,905,110 | 7,610,562 |
| AmpB_2 | 0.966007 | 0.964610 (6,370) | 0.933788 | 9,134,449 | 8,714,096 |
| AmpB_3 | 0.966353 | 0.964985 (6,372) | 0.934915 | 9,580,410 | 9,173,284 |
| pooled (38,862 values) | 0.964905 | 0.963385 (38,193) | 0.934965 | 54,832,659 | 52,253,712 |

Zero-count genes per sample: quant 119 to 152, align 261 to 289. Per gene across the 6 samples
(6,264 genes with nonzero variance in both modes): median Spearman 1 (6 values per gene: little
resolution).

DESeq2 (`~ condition`, Wald, contrast `condition,AmpB,Control`, no shrinkage), `padj < 0.05`,
NA padj excluded (`results/de_overlap.tsv`):

| measure | value |
|---|---|
| genes in each result | 6,477 (quant), 6,477 (align) |
| padj NA | 335 (quant), 427 (align) |
| padj < 0.05, quant | 2,042 |
| padj < 0.05, align | 2,019 |
| in both | 1,956 |
| quant only | 86 |
| align only | 63 |
| Jaccard (both / either) | 0.929216 |
| in both with the same log2FC sign | 1,954 of 1,956 |
| Spearman of log2FC, genes with a log2FC in both | 0.926467 |

No criterion was set for these numbers. I do not call them good or bad. There is no ground truth:
the GEO record lists no supplementary count table, and the authors counted CDS with HTSeq while
qRNAmap counts exons with featureCounts or kallisto. I did not compare with the authors' results.

### Bundling (measured)

The bundled settings ran: 6 samples in 1 quant task (23:29) and 3 samples per align task (20:51 and
22:13), all with `--threads 1`. The align task MaxRSS with 3 samples (13,070,508 and 12,895,776 KB)
is about the same as for 1 sample in E3 (13,073,276 KB), not 3 times larger. For a task of 1 to 1.5
hours at the rates measured here (quant 22.7 to 26.8 s, align 35.2 to 52.1 s per million pairs),
a task holds about 134 to 238 million pairs in quant mode and 69 to 153 million pairs in align
mode. These are arithmetic from this run, not a tested setting.

## E3. D2: why align tasks show a high MaxRSS

Sample: Control_2 (SRR3396386, 10,302,142 pairs, the largest of the six). All three jobs ran on node
r46 with `--mem 32G`, 1 CPU requested (2 allocated), at the same time as the E2 align tasks.

1. Manual steps (job 32308382, `run_amb.sh e3-manual Control_2`): the same commands as the CLI
   (`crates/qrnamap-core/src/toolcall.rs` `hisat2_align`, `feature_counts`; staging as in
   `crates/qrnamap-cli/src/pipeline.rs` `stage_reads`), with the index, splice sites and SAF of
   `$DATA/run_align`: `zcat` of each mate to a plain file in `$SCRATCH`;
   `hisat2-align-s -p 1 -x INDEX --known-splicesite-infile SS --no-temp-splicesite --max-intronlen
   5000 -1 R1.fq -2 R2.fq -S X.sam`; `featureCounts -T 1 -F SAF -a genes.saf -p --countReadPairs
   -s 0 -o featurecounts.tsv X.sam`. Each under `/usr/bin/time -v`. The featureCounts summary
   equals that of the CLI run 32308383 (only the SAM path in the header differs).
2. `qrnamap align --task-id 0 --samples-per-task 1` on Control_2 with the index copied from
   `$DATA/run_align` (`run_amb.sh e3-align Control_2 scratch|gpfs`): job 32308383 with
   `--tmpdir $SCRATCH/e3tmp`, job 32308384 with `--tmpdir $DATA/e3/tmp_gpfs_JOBID` (GPFS).
3. Cluster settings: `scontrol show config` (saved in `$DATA/e3/scontrol_show_config.txt`).

Measured (`results/e3_memory.tsv`):

| measure | value |
|---|---|
| `$SCRATCH` on r46 (`df -T`, `findmnt -T`) | xfs on `/dev/nvme0n1p6`, mounted at `/scratch` |
| staged plain FASTQ, mate 1 and mate 2 | 2,901,968,200 bytes each |
| SAM | 7,508,765,875 bytes |
| FASTQ + SAM | 13,312,702,275 bytes |
| `hisat2-align-s` Maximum resident set size | 41,208 KB (elapsed 5:39.74) |
| `featureCounts` Maximum resident set size | 36,144 KB (elapsed 0:17.80) |
| `zcat` Maximum resident set size | 3,340 KB and 3,260 KB |
| step cgroup `memory.stat` at job end: `total_cache` | 13,313,040,384 bytes |
| step cgroup `memory.stat` at job end: `total_rss` | 954,368 bytes |
| step cgroup `memory.max_usage_in_bytes` | 13,386,948,608 bytes |
| sacct MaxRSS, job 32308382 (manual) | 13,073,192 KB (x 1024 = 13,386,948,608 bytes) |
| sacct MaxRSS, job 32308383 (`qrnamap align`, `--tmpdir` on `$SCRATCH`) | 13,073,276 KB, elapsed 6:28 |
| sacct MaxRSS, job 32308384 (`qrnamap align`, `--tmpdir` on GPFS) | 39,392 KB, elapsed 7:18 |
| gene counts of 32308383 and 32308384 | identical (`cmp` of the decompressed files) |

Cluster accounting: `JobAcctGatherType = jobacct_gather/cgroup`, `JobAcctGatherFrequency = 30`,
`ProctrackType = proctrack/cgroup`, `ConstrainRAMSpace = yes`, `AllowedRAMSpace = 100.0%`. The job
memory cgroup is cgroup v1 (`/proc/self/cgroup` line `8:memory:/slurm_r46/uid_1181/job_ID/step_batch/task_0`;
the SLURM plugin directory has `cgroup_v1.so` and no `cgroup_v2.so`). `strings` of `cgroup_v1.so`
lists `memory.max_usage_in_bytes`, `memory.stat` and `total_rss`. `man sacct` says only "MaxRSS:
Maximum resident set size of all tasks in job."

What the measurements show:

- sacct MaxRSS of job 32308382 equals the step cgroup's `memory.max_usage_in_bytes` exactly
  (13,073,192 KB x 1024 = 13,386,948,608 bytes). On this cluster the sacct MaxRSS value of a job
  matched the cgroup v1 peak usage, which counts page cache. I observed this for one job.
- The cgroup's `total_cache` at the end (13,313,040,384 bytes) is within 338,109 bytes of the
  size of the staged FASTQ files plus the SAM (13,312,702,275 bytes). `total_rss` was 954,368 bytes.
- The largest process RSS was 41,208 KB (HISAT2). featureCounts used 36,144 KB.
- The same `qrnamap align` with its temporary files on GPFS had a MaxRSS of 39,392 KB, against
  13,073,276 KB with the files on `$SCRATCH` (xfs). Same sample, node, binary and output.
- 3 samples per align task gave the same MaxRSS as 1 sample (E2), which fits files that are deleted
  after each sample.

Hypotheses:

- A (accounting counts the page cache of the staged FASTQ and SAM in `$SCRATCH`): supported by every
  measurement above. The high MaxRSS is page cache of the files written to local scratch, not
  process memory.
- B (`$SCRATCH` is memory-backed, tmpfs): contradicted for nodes r46 and x01 (xfs on NVMe). The other nodes
  were not checked; the yeast align tasks ran on c07, which was not checked.
- C (one process holds the data): contradicted. HISAT2 41,208 KB, featureCounts 36,144 KB.

Bytes per pair for Control_2: 13,386,948,608 / 10,302,142 = 1,299 bytes per pair (yeast: 1,329 to
1,338).

Memory limit below the peak: job 32309177 ran the same `run_amb.sh e3-align Control_2 scratch` with
`--mem 4G`. It ran on node x01 (`$SCRATCH` there: xfs on `/dev/nvme0n1p6`). It COMPLETED with exit
0:0 in 10:32, MaxRSS 4,192,060 KB (the 4G limit is 4,194,304 KB). Its gene counts are identical to
those of job 32308383 (`cmp`). So the 13 GB page-cache peak did not need a 13 GB `--mem` in this one
run: the cgroup kept usage at the limit and the job finished. The elapsed time (10:32 against 6:28)
came from another node, so the run does not show whether the limit slowed the job.

## Defects found

D1 (from the yeast run): the generated scripts no longer depend on the PATH of the submit shell.
Verified on real jobs in E2 (baked PATH worked). The cancel of dependent jobs after a failed prepare
was observed (E1).

D2 (from the yeast run): explained by E3: the align MaxRSS is page cache of the scratch files
counted by the cgroup accounting, not process memory. Not changed in code. One run with `--mem 4G`
completed with identical counts (job 32309177). Open: a limit below 4 GB, other nodes, and the
effect of the limit on run time.

D3 (new). **The generated SLURM scripts need `$SCRATCH`, which on this cluster comes from the
submitting shell's module environment.** The jobs inherit the environment of `sbatch` (default
`--export=ALL`). The SLURM task prolog
(`/opt/linux/rhel/8.x/x86_64/pkgs/slurm/24.11.1/scripts/task-prolog.sh`) runs `modulecmd bash
reload`, and the module `workspace/scratch` sets `SCRATCH=/scratch/$USER/$SLURM_JOB_ID`. When the
submit shell has no loaded modules, the job gets no `SCRATCH` and every generated job fails at its
first line `: "${SCRATCH:?...}"`.
- Failing command: `env -i HOME="$HOME" USER="$USER" PATH="<slurm bin>:/usr/bin:/bin" bash $DATA/e1/run_quant/slurm/submit.sh`
- Message (job 32308072): `/var/spool/slurmd/job32308072/slurm_script: line 17: SCRATCH: SCRATCH is not set: this job needs node-local scratch`
- Minimal reproduction without qrnamap: from the same `env -i` shell,
  `sbatch -p short -t 00:02:00 --mem 100M --wrap 'echo SCRATCH=${SCRATCH:-UNSET}'` prints
  `SCRATCH=UNSET` (job 32308141 had `LOADEDMODULES=workspace/scratch` set as well; the module reload
  failed because the module needs `HPCC_MODULES`). From a login shell (`bash -l`) the jobs get
  `SCRATCH` (E1b, E2).
- Effect: the job fails fast and with a clear message; with `--kill-on-invalid-dep=yes` the
  dependent jobs are cancelled (E1a). No wrong output. The R28 note in the scripts ("the jobs do not
  depend on the PATH of the shell that runs submit.sh") is true for PATH, but the jobs still depend
  on the submit environment for `$SCRATCH`. Not fixed (no CLI change in this task).

## Not measured

- No ground truth: the numbers show agreement between the two modes, not accuracy. No comparison
  with the authors' HTSeq counts (none published on the GEO sample pages) or their DE results.
- Local mode (`qrnamap run` without `--slurm`) on this data.
- `--threads` above 1; `--keep-sam` in a CLI run; `--shrink`.
- `$SCRATCH` file system type on nodes other than r46 and x01 (the yeast align node c07 was not checked).
- Process RSS of the `qrnamap` process itself during staging (only the tools were wrapped in
  `/usr/bin/time -v`).
- Why the later samples of each align task took longer than the first (r46 was shared with my
  other jobs from 09:52 to 10:00).
- The per-sample seconds are derived from file mtimes, not from timers.
- strand-check used only the first 100,000 pairs per sample (default).
- Whether the TruSeq RNA v2 kit is unstranded is not in the record; strand-check measured it.
- Which sacct field the accounting plugin uses was inferred from one exact equality, not from the
  SLURM source.
