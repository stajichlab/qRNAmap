# qRNAmap

A Rust command-line tool for bulk RNA-seq quantification: genome + GFF3 + sample sheet + FASTQ to
gene counts (kallisto, or HISAT2 + featureCounts) and DESeq2. Licence: GPL-3.0-only.

## Install and build

The tools (kallisto, HISAT2, featureCounts, gffread, R with DESeq2, Rust) come from the pixi
environment in `pixi.toml`.

```bash
pixi install                            # not tested on the UCR login node
cargo build --release -p qrnamap-cli    # with the pixi env on PATH, see below
# the binary is target/release/qrnamap
```

On the UCR HPCC login node, `pixi run ...` hangs (cause not found). Put the environment on PATH
yourself instead:

```bash
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
cargo build --release -p qrnamap-cli
QRNAMAP_REQUIRE_TOOLS=1 cargo test --workspace
```

## Sample sheet

A tab-separated file with a header row.

| column | content |
|---|---|
| `sample` | required. The sample name. It starts with a letter and uses only `[A-Za-z0-9_.-]`. |
| `fastq1` | required. One FASTQ path, or a comma-separated list of lanes. Plain, `.gz` or `.zst`. |
| `fastq2` | optional. Mate 2, with the same number of lanes as `fastq1`. Empty means single-end. |
| any other column | a covariate. A factor by default. A header `name:numeric` makes a numeric column `name`. |

- Covariate column names start with a letter and use only `[A-Za-z0-9_.]`. R reserved words are refused.
- Factor values use only `[A-Za-z0-9_.-]`, with no space. Write `wild_type`, not `wild type`.
- Blank lines and lines that start with `#` are ignored.
- Single-end samples need `--fragment-length MEAN,SD` in quant mode.

Example:

```
sample	fastq1	fastq2	condition
WT_1	/bigdata/x/WT_1_L1_R1.fq.gz,/bigdata/x/WT_1_L2_R1.fq.gz	/bigdata/x/WT_1_L1_R2.fq.gz,/bigdata/x/WT_1_L2_R2.fq.gz	wild_type
KO_1	/bigdata/x/KO_1_R1.fq.gz	/bigdata/x/KO_1_R2.fq.gz	ko
```

## Commands

```bash
# Genome FASTA and GFF3 (NCBI from the built-in catalog, or FungiDB with an API key)
qrnamap fetch --source ncbi --id GCF_000146045.2 --out genome/
qrnamap fetch --source fungidb --id "Aspergillus fumigatus Af293" --out genome/ --eupathdb-key-file ~/.eupathdb_key

# Transcript FASTA and tx2gene from a genome and GFF3
qrnamap transcripts --reference genome.fna.gz --gff genes.gff.gz --out tx.fa.zst --tx2gene tx2gene.tsv.zst \
  --genes-table genes.tsv.zst

# kallisto per sample (writes OUT/<sample>.genecounts.tsv.zst)
qrnamap quant --samples samples.tsv --transcripts tx.fa.zst --tx2gene tx2gene.tsv.zst --out run/ --strandedness reverse

# HISAT2 + featureCounts per sample
qrnamap align --samples samples.tsv --reference genome.fna.gz --gff genes.gff.gz --out run/ --strandedness reverse

# Gene-by-sample matrix from the per-sample count files
qrnamap matrix --out run/ --samples samples.tsv --genes-from tx2gene tx2gene.tsv.zst

# DESeq2
qrnamap de --counts run/counts.tsv.zst --samples samples.tsv --contrast condition,ko,wild_type --out run/deseq2

# Strandedness from the first 100,000 pairs of each sample
qrnamap strand-check --samples samples.tsv --reference genome.fna.gz --gff genes.gff.gz --out strand/

# All steps on this machine
qrnamap run --samples samples.tsv --reference genome.fna.gz --gff genes.gff.gz --mode quant \
  --out run/ --strandedness reverse --contrast condition,ko,wild_type

# SLURM scripts instead of running (the run directory must be on shared storage).
# The jobs put temporary files in node-local $SCRATCH: run submit.sh from a login shell.
qrnamap run --samples samples.tsv --reference genome.fna.gz --gff genes.gff.gz --mode align \
  --out /bigdata/lab/run/ --strandedness reverse --contrast condition,ko,wild_type \
  --slurm --tool-path "$PWD/.pixi/envs/default/bin" --time 02:00:00
bash /bigdata/lab/run/slurm/submit.sh
```

The temporary directory is `--tmpdir`, then `$SCRATCH`, then `$TMPDIR`. A run directory under
`$SCRATCH` or `/scratch/` gives a warning: SLURM array tasks on other nodes cannot see it.

### SLURM and the temporary directory

- Default (no `--tmpdir`): each job puts its temporary files in node-local `$SCRATCH`. Submit
  `submit.sh` from a login shell. A job that was submitted from a bare environment has
  no `$SCRATCH` and stops at its first check with this message: `SCRATCH is not set: this job needs
  node-local scratch. Submit from a login shell where $SCRATCH is set, or regenerate the scripts with
  --tmpdir DIR.` Evidence: `validation/amb/README.md`, D3.
- Fallback for a submit from a bare environment: `run --slurm --tmpdir DIR` writes `DIR` into the
  scripts. The jobs then do not need `$SCRATCH`. A job creates `DIR` (and a `qrnamap-<pid>`
  subdirectory in it) when the parent of `DIR` is writable. Use storage that every compute node can
  reach, for example `/bigdata`. A `--tmpdir` under `$SCRATCH` or `/scratch` gives a warning, because
  such a path is job-specific. It is not an error.
- The cost of `--tmpdir` on shared storage: the staged FASTQ files and the SAM are written there. For
  one AMB align sample of 10,302,142 pairs, the staged FASTQ plus the SAM were 13,312,702,275 bytes
  (`validation/amb/README.md`, E3). qrnamap removes these files when a sample ends, also after an
  error. A job that SLURM stops (time limit, `scancel`, out of memory) cannot run this cleanup, so its
  files stay in `DIR`; remove them by hand. This follows from how the cleanup works (it runs in the
  program); it was not tested.
- Process IDs can repeat on different nodes, so two jobs on a shared `--tmpdir` can get the same
  `qrnamap-<pid>` subdirectory. The per-sample subdirectories in it keep the files of different
  samples apart. No collision was found.
- `--mem-gb` default: 8 (GB). The value is an unmeasured placeholder. Reason, from `validation/amb/README.md`
  (E3, one align sample of 10,302,142 pairs): `sacct` MaxRSS of an align job was 13,073,276 KB with the
  temporary files on `$SCRATCH` and 39,392 KB with them on GPFS (`/bigdata`). The 13 GB was, in the one job
  examined, the page cache of the staged files in `$SCRATCH` that the job accounting counts. The tools themselves were small: HISAT2 used 41,208 KB
  and featureCounts 36,144 KB resident. A job with `--mem 4G` finished with identical gene counts
  (`sacct` MaxRSS 4,192,060 KB on node x01 and 4,191,164 KB on node c17). That is two runs of one sample on
  two nodes. No job was run at 8 GB (the tested limits were 32 GB and 4 GB), and the two 4 GB runs took
  10:32 and 19:26 against 6:28 at 32 GB on different nodes, so the effect of the limit on run time is not
  attributed and stays open (`validation/amb/README.md`, E3). The yeast align tasks peaked at 6.7 to 8.0 GB
  MaxRSS (`validation/README.md`), at or near 8 GB. It was NOT measured for genomes near 50 Mb or for 20
  million pairs per sample: raise `--mem-gb` for those. The `--time` and `--cpus` defaults are still
  unmeasured placeholders. Measured AMB rates with `--threads 1`: quant 22.7 to 26.8 s and align 35.2 to
  52.1 s per million read pairs (`validation/amb/README.md`, "Bundling").

## Gene table and mapping summaries

These two files are information only. The count matrix `counts.tsv.zst` is unchanged, and nothing
is excluded from DESeq2.

### `genes.tsv.zst`

One row per gene, tab-separated, with a header. Columns: `gene contig strand biotype class mito n_transcripts`.

- `biotype`: the `gene_biotype` attribute of the gene in the GFF3. If the gene has none: `mRNA` when any
  transcript of the gene is `mRNA`; else the type of the first transcript; else `unknown`.
- `class`: one of `protein_coding`, `rRNA`, `tRNA`, `other_ncRNA`, `pseudogene`, `other`, derived from
  `biotype` (case-insensitive) by these rules in this order: `protein_coding` or `mRNA` gives
  `protein_coding`; `rRNA` gives `rRNA`; `tRNA` gives `tRNA`; a biotype that contains `pseudogene` or
  `pseudogenic` gives `pseudogene`; any other biotype that ends with `rna` gives `other_ncRNA`; everything
  else gives `other`.
- `mito`: `yes` or `no`. A gene is `yes` when its contig has a `region` feature with the attribute
  `genome=mitochondrion` or `genome=mitochondrial` in the GFF3 (the match is case-sensitive), or when the
  contig is named with `--mito-contig NAME` (repeatable). The two sources are added together. The
  `--mito-contig` match is exact and case-sensitive. A name that matches no contig of the GFF3 gives a
  note on stderr, not an error.
- `n_transcripts`: the number of transcripts of the gene.

Which subcommands write it:

- `qrnamap transcripts --genes-table PATH` (the compression follows the file suffix).
- `qrnamap quant --gff GFF`: writes `OUT/genes.tsv.zst`. Plain `quant` without `--gff` writes no table
  (stderr says so). `--mito-contig` without `--gff` on plain `quant` is an error. Every gene of the
  `--tx2gene` table must be a gene of the `--gff` GFF3, else `quant` stops with an error that names the
  missing genes, before kallisto runs. A GFF3 with more genes is accepted.
- `quant --task-id` and `align --task-id` do not write the table (the prepare step does); `--mito-contig`
  there gives a note that it has no effect. `transcripts --mito-contig` without `--genes-table` gives the
  same note.
- `qrnamap align` and `qrnamap run`: write `OUT/genes.tsv.zst` from the required `--gff`.
- `qrnamap strand-check` does not write or change `genes.tsv.zst`, also when its `--out` is a run directory.

The table is written after the index is built, so a failed prepare step leaves no table. A rerun
removes an old table first.

### `<sample>.summary.tsv.zst` and `mapping_summary.tsv.zst`

`quant` and `align` write `OUT/<sample>.summary.tsv.zst` as the last step for each sample. It has a
header and one row. Columns: `sample mode reads_processed fraction_mapped fraction_assigned strandedness qrnamap_version`.
`mode` is `quant` or `align`.

`matrix` (and `run`, and `run --step final`) read these files in sample-sheet order and write
`OUT/mapping_summary.tsv.zst`: the same columns plus `strand_check_suggestion`.

- `strand_check_suggestion` is the `suggested` cell of `OUT/strand_check.tsv` for the sample. The file
  is read only when it is in the same directory as `matrix --out`. The value is copied as it is. It is
  not checked against the strandedness of the run, and it carries no stamp. Without the file, the
  column is empty.
- Floats are printed in the shortest form that reads back to the same value. `1.0` prints as `1`. A
  column in which every value is 1 is therefore read by pandas as integer.
- `matrix` stops with a named error when the summary of a listed sample is missing or belongs to
  another sample. Run directories made before this version have no summary files: rerun each
  sample (`quant` or `align`) before `matrix` or `run --step final`. Before any check, `matrix` (and
  `run --step final`) removes the old `counts.tsv.zst` and `mapping_summary.tsv.zst`, so a failed `matrix`
  leaves neither. It does not remove old `deseq2.*` files; `de` removes those itself.
- `strand-check` removes an old `strand_check.tsv` in its `--out` after its input checks and before it
  runs any tool.

### Rounding of kallisto counts

The count of a gene from kallisto is the sum of `est_counts` over the transcripts of the gene. The sum is
rounded once per gene, half away from zero. It is not rounded per transcript.

## Output files (in the run directory)

| file | content |
|---|---|
| `<sample>.genecounts.tsv.zst` | gene counts of one sample (`gene`, sample name) |
| `<sample>/abundance.tsv`, `<sample>/run_info.json` | kallisto output (quant mode) |
| `<sample>/featurecounts.summary` | featureCounts summary (align mode) |
| `<sample>.summary.tsv.zst` | mapping summary of one sample (written by `quant` and `align`) |
| `counts.tsv.zst` | gene-by-sample matrix |
| `mapping_summary.tsv.zst` | one row per sample, from the per-sample summaries |
| `genes.tsv.zst` | gene table: class, biotype and mitochondrial flag per gene |
| `deseq2.results.csv`, `deseq2.normalized.csv`, `deseq2.vst.csv`, `deseq2.versions.txt` | DESeq2 output, with `--contrast` |
| `deseq2.results.shrunk.csv` | with `--shrink` |
| `transcripts.fa.zst`, `tx2gene.tsv.zst`, `transcripts.kidx` | reference files (quant mode) |
| `splice_sites.txt`, `genes.saf`, `index/` | reference files (align mode) |
| `.prepared` | marker of the prepare step; array tasks and the final step check it |
| `slurm/` | the SLURM scripts and logs, with `--slurm` |
| `strand_check.tsv` | from `strand-check`, in its `--out` |

## Defaults that are not measured

These values are chosen, not measured or calibrated. Check them for your data.

- HISAT2 `--max-intronlen` 5000 (`--max-intron`). Not verified for fungi.
- SLURM `--time 02:00:00`, `--mem-gb 8`, `--cpus` (default: the `--threads` value). See "SLURM and the temporary directory" for the `--mem-gb` evidence.
- strand-check threshold 0.8: a strand is suggested when its fraction is above 0.8.
- `LOW_ASSIGNED_FRACTION` 0.5: a sample with an assigned (or pseudoaligned) fraction below 0.5 gives a warning.

Align memory, measured on the yeast validation (`validation/README.md`): 6.7 to 8.0 GB MaxRSS per
array task at 5.2 to 6.1 million pairs with `--threads 1`. This is consistent with page cache of the
staged scratch files, as measured in the AMB experiments (`validation/amb/README.md`, E3; one job showed the equality of
`sacct` MaxRSS and the cgroup peak that counts page cache).

## Not in this version

- No browser page and no WASM build of the tool. Only the three pure crates are compiled for
  `wasm32-unknown-unknown` as a check.
- No Salmon.
- No plots.
- `--shrink apeglm` and `--shrink ashr` need the R packages apeglm and ashr. `pixi.toml` does not
  install them, so these options have never run to success here. Only `--shrink normal` was run.
- The SLURM `--time` is not computed from `reads_per_sample`.

## Validation

Two validation records are in the repository. Both ran on the UCR HPCC cluster.

- `validation/README.md`: yeast, PRJNA674653 (6 paired-end S288c samples, 5.2 to 6.1 million pairs each). Headline numbers from that file:
  - Transcript FASTA identical to gffread for all 6,478 transcripts.
  - strand-check suggested `reverse` for all 6 samples (antisense fraction 0.9855 to 0.9880).
  - Quant: 93.10 to 94.60% pseudoaligned. Align: 96.65 to 97.53% HISAT2 overall alignment, assigned fraction 0.8875 to 0.8979.
  - Spearman correlation of gene counts between the two modes: 0.971456 to 0.973390 per sample.
  - DESeq2 genes with padj < 0.05: 2,648 (quant), 2,632 (align), 2,571 in both (Jaccard 0.949059).
- `validation/amb/README.md`: AMB, GEO GSE80357 (6 paired-end S. cerevisiae samples, 8,405,534 to 10,302,142 pairs). Headline numbers from that file:
  - strand-check suggested `unstranded` for all 6 samples (antisense fraction 0.4951 to 0.4988).
  - Quant: 93.6 to 94.7% pseudoaligned. Align: 95.49 to 97.06% HISAT2 overall alignment, assigned fraction 0.8802 to 0.9087.
  - Spearman correlation of gene counts between the two modes, pooled over all 6 samples: 0.964905.
  - DESeq2 genes with padj < 0.05: 2,042 (quant), 2,019 (align), 1,956 in both (Jaccard 0.929216).
  - The `--kill-on-invalid-dep` flag cancelled the dependent jobs 12 s after a failed prepare job, in 2 of 2 attempts.
- No ground-truth comparison was made. For GSE80357 the GEO series page lists the authors' counts file and
  fold-change table. They were not downloaded and not compared. The numbers show agreement between the two modes, not accuracy.

## Layout

- `qseq-io`, `genome-catalog`, `qrnamap-core`: pure libraries (no file or network I/O). They must build for `wasm32-unknown-unknown`.
- `genome-fetch`: library for downloads.
- `qrnamap-cli`: the `qrnamap` binary.

## Development

- `pixi run test` runs `cargo test --workspace` with `QRNAMAP_REQUIRE_TOOLS=1`.
- `pixi run wasm-check` builds the three pure crates for `wasm32-unknown-unknown`.

## Tool versions

All packages solve in one pixi environment (`default`). Solved on 2026-11-19 (linux-64):

| Package | Version | Channel |
|---|---|---|
| kallisto | 0.52.0 | bioconda |
| hisat2 | 2.2.2 | bioconda |
| subread (featureCounts) | 2.1.1 | bioconda |
| gffread | 0.12.9 | bioconda |
| r-base | 4.5.3 | conda-forge |
| bioconductor-deseq2 | 1.50.2 | bioconda |
| rust | 1.100.0 (constraint `>=1.96,<2`) | conda-forge |

## Measured

Pipeline recovery on simulated reads (Task 11, 2026-11-19). Tests: `crates/qrnamap-cli/tests/pipeline_quant.rs`, `pipeline_align.rs`.

- Data: 10,000 error-free read pairs from the 6 transcripts of `tests/fixtures/mini` (5 genes, 1,000 to 3,000 pairs per gene). The transcripts are 70 to 120 bases long, shorter than the 200-base fragment, so each fragment is the whole transcript. Reads are 75 bases (70 for the 70-base transcript).
- kallisto 0.52.0 (`-t 1`): every gene count equals the truth (maximum relative error 0). Test tolerance: 10%.
- HISAT2 2.2.2 + featureCounts 2.1.1 (`-p 1`, `-T 1`): every gene count equals the truth (maximum relative error 0); all 10,000 pairs mapped and assigned. Test tolerance: 5%.
- Wrong strandedness: a reverse library run as forward gives featureCounts assigned fraction 0.0 and kallisto 0 pseudoaligned reads (kallisto then exits with status 1). An unstranded library run as forward in kallisto gives a pseudoaligned fraction of 0.496.
- These results come from very short, error-free transcripts. They do not predict accuracy on real data.

Native HISAT2 runs `hisat2-align-s`, not the `hisat2` wrapper: the Perl wrapper passes its arguments through a shell, which split a path with a space and expanded `$HOME` in the tests.

Larger synthetic fixture (`crates/qrnamap-cli/tests/common/large.rs`, test `tests/large.rs`, 2026-11-19): 2 contigs of 20 kb, 12 genes (2 multi-isoform, 1 rRNA, 1 CDS-only), 15 transcripts of 661 to 2,034 bases, 40,000 error-free pairs (read length 100, fragment 200 to 300). The reads come from the generator's own transcripts, not from our builder. Measured maximum relative error per gene (genes with at least 200 pairs): kallisto 0 (unstranded, forward and reverse libraries with the matching flag); HISAT2 + featureCounts 0.0005 unstranded and 0.0018 stranded. The lost pairs are g07 pairs that map to 2 places (`Unassigned_MultiMapping`). A reverse library run as forward: kallisto pseudoaligns no read (exit status 1); featureCounts assigned fraction 0.0.
