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
qrnamap transcripts --reference genome.fna.gz --gff genes.gff.gz --out tx.fa.zst --tx2gene tx2gene.tsv.zst

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

# SLURM scripts instead of running (the run directory must be on shared storage)
qrnamap run --samples samples.tsv --reference genome.fna.gz --gff genes.gff.gz --mode align \
  --out /bigdata/lab/run/ --strandedness reverse --contrast condition,ko,wild_type \
  --slurm --tool-path "$PWD/.pixi/envs/default/bin" --mem-gb 16 --time 02:00:00
bash /bigdata/lab/run/slurm/submit.sh
```

The temporary directory is `--tmpdir`, then `$SCRATCH`, then `$TMPDIR`. A run directory under
`$SCRATCH` or `/scratch/` gives a warning: SLURM array tasks on other nodes cannot see it.

## Output files (in the run directory)

| file | content |
|---|---|
| `<sample>.genecounts.tsv.zst` | gene counts of one sample (`gene`, sample name) |
| `<sample>/abundance.tsv`, `<sample>/run_info.json` | kallisto output (quant mode) |
| `<sample>/featurecounts.summary` | featureCounts summary (align mode) |
| `counts.tsv.zst` | gene-by-sample matrix |
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
- SLURM `--time 02:00:00`, `--mem-gb 16`, `--cpus` (default: the `--threads` value).
- strand-check threshold 0.8: a strand is suggested when its fraction is above 0.8.
- `LOW_ASSIGNED_FRACTION` 0.5: a sample with an assigned (or pseudoaligned) fraction below 0.5 gives a warning.

Align memory, measured on the yeast validation (`validation/README.md`): 6.7 to 8.0 GB MaxRSS per
array task at 5.2 to 6.1 million pairs with `--threads 1`. The cause was not found. The 16 GB
default `--mem-gb` is not measured for larger inputs; set `--mem-gb` yourself.

## Not in this version

- No browser page and no WASM build of the tool. Only the three pure crates are compiled for
  `wasm32-unknown-unknown` as a check.
- No Salmon.
- No plots.
- `--shrink apeglm` and `--shrink ashr` need the R packages apeglm and ashr. `pixi.toml` does not
  install them, so these options have never run to success here. Only `--shrink normal` was run.

## Validation on public yeast RNA-seq

See `validation/README.md` (6 paired-end S288c samples, 5.2 to 6.1 million pairs each, PRJNA674653).
Headline numbers from that file:

- Transcript FASTA identical to gffread for all 6,478 transcripts.
- strand-check suggested `reverse` for all 6 samples (antisense fraction 0.9855 to 0.9880).
- Quant: 93.10 to 94.60% pseudoaligned. Align: 96.65 to 97.53% HISAT2 overall alignment, assigned fraction 0.8875 to 0.8979.
- Spearman correlation of gene counts between the two modes: 0.971456 to 0.973390 per sample.
- DESeq2 genes with padj < 0.05: 2,648 (quant), 2,632 (align), 2,571 in both (Jaccard 0.949059).
- No ground truth exists for this data: these numbers show agreement, not accuracy.

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
