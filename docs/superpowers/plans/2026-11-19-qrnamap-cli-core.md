# qRNAmap CLI and Core Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A working native `qrnamap` CLI that goes from genome + GFF3 + a sample sheet + FASTQ to a gene count matrix and a DESeq2 result, in fast mode (kallisto) and genome mode (HISAT2 + featureCounts).

**Architecture:** A Cargo workspace of small crates. Pure crates (`qseq-io`, `genome-catalog`, `qrnamap-core`) do no I/O and compile for `wasm32`. `qrnamap-cli` and `genome-fetch` are native. Native tools run through a serialisable `ToolCall`.

**Tech Stack:** Rust (pixi `rust >=1.96,<2`), pixi, kallisto 0.52.0, HISAT2 2.2.2, Subread/featureCounts 2.1.1, gffread, R + DESeq2, pytest for the catalog generator.

**Spec:** `docs/superpowers/specs/2026-11-19-qrnamap-design.md` (approved 2026-11-19). Section numbers below refer to it.

**Scope of this plan (CLI first, spec 11a).** This plan covers build steps 1 to 4 and the CLI side of step 7. It does not cover the browser. Spike 0 (all of 0a to 0f), the `qrnamap-wasm` crate, the web page, WASM lowering of `ToolCall`, chunked streaming and the Salmon backend belong to Plan 2 and Plan 3. They are written after Plan 1, using what Plan 1 and the spike measure.

## Global Constraints

- Licence: GPL-3.0-only; `LICENSE` is already in the repo (spec 11a).
- pixi manifests use `[workspace]`, never `[project]` (global rule).
- `qseq-io`, `genome-catalog`, `qrnamap-core` have no file or network I/O and must build with `cargo build --target wasm32-unknown-unknown` (spec 5).
- Tool versions: kallisto 0.52.0, HISAT2 2.2.2, featureCounts 2.1.1 (spec 3). DESeq2: bioconda has 1.50.2; start with it, print the version in every DE output, and swap to 1.52.0 from Bioconductor source later (spec 11c).
- Large text output is compressed: `zstd` for internal files, `gzip` for shared files. New readers accept plain, `.gz` and `.zst` (global rule, spec 11f).
- Job-local temporary files go in `${SCRATCH:?}`; never hardcode `/scratch/...`. Never use `BASH_SOURCE[0]` in scripts that run under SLURM. Index files go on `/bigdata` before an array starts.
- `reads_per_sample` is a parameter, default 20,000,000 pairs (spec 11a).
- Native parity runs use `kallisto -t 1` and `hisat2 -p 1`. Gene counts from kallisto: sum `est_counts` per gene, round to integer (spec 6).
- Default DE: `~ condition`, Wald test, no LFC shrinkage (spec 11b). Design strings use a strict grammar and are quoted in generated R (spec 11c).
- Strand `.` or `?` in GFF3 is an error, not plus. A `##FASTA` section must parse. All isoforms are kept (spec 11c).
- Featurecounts runs with `-F SAF` made from our gene model, with `-p --countReadPairs` for pairs. Both modes must give identical gene row names.
- Never report speed or accuracy numbers that were not measured.

## Review Focus

1. A GFF3 with a `##FASTA` section, CRLF line endings, or `.gz` compression: parses; never silently drops genes. (Task 3, Task 4)
2. Two genes with the same label, or a transcript on a contig absent from the FASTA: a named error, not a wrong matrix. (Task 3, Task 5)
3. A single-end sample with no fragment length mean and SD for kallisto: a named error before any tool runs. (Task 9)
4. A sample name or factor level with a space, quote or `$` that would reach the generated R script: rejected or quoted, never executed. (Task 8, Task 12)
5. Fewer than 2 samples in a tested group, or an all-zero sample, in DE: refuse with the reason. (Task 12)
6. A wrongly set strandedness (most reads assigned to the opposite strand): reported as a low assigned fraction per sample. (Task 11)

## File Structure

```
Cargo.toml                       workspace (members below)
pixi.toml                        [workspace]; tools and Rust toolchain
.github/workflows/ci.yml         cargo test + wasm32 build of pure crates
crates/qseq-io/                  FASTA/FASTQ readers; GFF3 gene model (all isoforms)
crates/genome-catalog/           catalog schema, URLs, size checks; generator script + pytest
crates/genome-fetch/             native downloads (NCBI, FungiDB with key)
crates/qrnamap-core/             transcripts, SAF, splice sites, sample sheet, design, ToolCall,
                                 count parsers, matrix, DE script, strand check
crates/qrnamap-cli/              binary `qrnamap`: IO, executors, subcommands
tests/fixtures/mini/             synthetic genome, GFF3, truth counts (made in Task 5)
```

---

### Task 1: Workspace, pixi environment, CI

**Files:**
- Create: `Cargo.toml`, `pixi.toml`, `.gitignore`, `.github/workflows/ci.yml`, `README.md`, and an empty `lib.rs` (or `main.rs`) with one test in each crate listed above.

**Interfaces:**
- Produces: workspace members `qseq-io`, `genome-catalog`, `genome-fetch`, `qrnamap-core`, `qrnamap-cli`; the pixi tasks `test` (`cargo test --workspace`) and `wasm-check` (builds the three pure crates for `wasm32-unknown-unknown`).

- [ ] **Step 1: Write a smoke test in each crate** named `crate_builds` (`assert!(true)`), and a pixi task `wasm-check` that builds `-p qseq-io -p genome-catalog -p qrnamap-core --target wasm32-unknown-unknown`.
- [ ] **Step 2: Write `pixi.toml`** with `[workspace]`, channels `conda-forge` and `bioconda`, platform `linux-64`, dependencies `rust >=1.96,<2`, `rust-std-wasm32-unknown-unknown >=1.96,<2`, `kallisto ==0.52.0`, `hisat2 ==2.2.2`, `subread ==2.1.1`, `gffread`, `zstd`, `python >=3.11`, `pytest`, `r-base`, `bioconductor-deseq2 ==1.50.2`. Run `pixi install` and record the solved versions of kallisto, hisat2, subread, gffread, r-base in `README.md` under "Tool versions".
- [ ] **Step 3: Run** `pixi run test` and `pixi run wasm-check`. Expected: both pass.
- [ ] **Step 4: Write `ci.yml`** running the same two tasks with `prefix-dev/setup-pixi` on push and pull request.
- [ ] **Step 5: Commit** `git commit -m "chore: workspace, pixi environment, CI"`; push; confirm the CI run is green with `gh run watch`.

### Task 2: `qseq-io` sequence readers

**Files:**
- Create: `crates/qseq-io/src/lib.rs`, `crates/qseq-io/src/fasta.rs`, `crates/qseq-io/src/fastq.rs`
- Test: `crates/qseq-io/tests/readers.rs`

**Interfaces:**
- Produces (copied from qmutfind `qmutfind-core/src/parse.rs`, same behaviour; record the qmutfind commit hash in the commit message):
  `pub fn read_fasta<R: BufRead>(reader: R) -> io::Result<Vec<(String, Vec<u8>)>>`
  `pub struct FastqReader<R: BufRead>` with `new(reader)` and an iterator of records
  `pub struct FastqChunker` with `new()`, `with_limit(usize)`, `push(&mut self, &[u8]) -> io::Result<Vec<(String, Vec<u8>)>>`, `finish(&mut self) -> io::Result<Vec<(String, Vec<u8>)>>` (kept now so the browser plan can reuse it).
- Compression detection is not in this crate (spec 11c).

- [ ] **Step 1: Copy qmutfind's `qmutfind-core/tests/chunker.rs` and the FASTA/FASTQ cases** into `tests/readers.rs`, renamed to this crate's paths. Add `fasta_multiline_and_lowercase_kept` (a two-line sequence joins; lowercase stays lowercase).
- [ ] **Step 2: Run** `pixi run cargo test -p qseq-io`. Expected: FAIL (functions missing).
- [ ] **Step 3: Copy `parse.rs` into `fasta.rs` and `fastq.rs`**, no behaviour change; re-export from `lib.rs`.
- [ ] **Step 4: Run** `pixi run cargo test -p qseq-io` and `pixi run wasm-check`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(qseq-io): FASTA and FASTQ readers extracted from qmutfind@<hash>`.

### Task 3: `qseq-io` gene model (all isoforms)

**Files:**
- Create: `crates/qseq-io/src/gff.rs`
- Test: `crates/qseq-io/tests/gff.rs`

**Interfaces:**
- Produces:
  ```rust
  pub enum Strand { Plus, Minus }
  pub struct Transcript { pub id: String, pub kind: String,            // mRNA, rRNA, tRNA, ncRNA, ...
      pub exons: Vec<(u64,u64)>, pub cds: Vec<(u64,u64)>,              // 1-based inclusive, sorted ascending
      pub utr5: Vec<(u64,u64)>, pub utr3: Vec<(u64,u64)> }             // explicit UTR features only
  pub struct Gene { pub label: String, pub contig: String, pub strand: Strand, pub transcripts: Vec<Transcript> }
  pub struct GeneModel { pub genes: Vec<Gene> }
  pub enum GffError { BadLine { line: usize }, BadStrand { line: usize, value: String },
      DuplicateGeneLabel(String), Io(io::Error) }
  pub fn parse_gff3<R: BufRead>(reader: R) -> Result<GeneModel, GffError>
  ```
- Gene label rule: `locus_tag`, else `Name`, else `ID` (as qmutfind `annotate.rs` lines 236-241). A transcript's gene is its `Parent`; a transcript with no gene parent becomes its own gene.

- [ ] **Step 1: Write tests** in `tests/gff.rs`, each with a small inline GFF3 string:
  `keeps_all_isoforms` (gene with two mRNAs → 2 transcripts); `stops_at_fasta_section` (lines after `##FASTA` ignored, parse succeeds); `crlf_lines_parse`; `strand_dot_is_error` (→ `GffError::BadStrand{line, value:"."}`); `short_line_is_error` (→ `BadLine`); `rrna_trna_kept_with_kind` (`kind == "rRNA"`); `utr_features_collected`; `duplicate_gene_label_is_error`; `label_prefers_locus_tag_over_name_over_id`.
- [ ] **Step 2: Run** `pixi run cargo test -p qseq-io --test gff`. Expected: FAIL.
- [ ] **Step 3: Implement `parse_gff3`** as a two-pass parse (collect features by `ID`, then attach children by `Parent`). New code; do not copy qmutfind's canonical-transcript logic.
- [ ] **Step 4: Run** the tests and `pixi run wasm-check`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(qseq-io): GFF3 gene model with all isoforms`.

### Task 4: Compressed input in the CLI

**Files:**
- Create: `crates/qrnamap-cli/src/io.rs`
- Test: `crates/qrnamap-cli/tests/io.rs`

**Interfaces:**
- Produces: `pub fn open_input(path: &Path) -> Result<Box<dyn BufRead>, InputError>` (plain, `.gz` multi-member via `flate2::bufread::MultiGzDecoder`, `.zst` via `zstd::Decoder`); `pub enum InputError { Io(io::Error), Truncated { path: PathBuf } }`. Detection is by file extension and then by magic bytes (`1f 8b`, `28 b5 2f fd`); a mismatch uses the magic bytes.
- Consumes: nothing.

- [ ] **Step 1: Write tests** (`tempfile` files): `reads_plain`, `reads_gzip_single_member`, `reads_gzip_multi_member` (two members concatenated), `reads_zstd`, `extension_wrong_magic_wins` (a `.txt` that holds gzip), `truncated_gzip_is_named_error` (→ `InputError::Truncated`).
- [ ] **Step 2: Run** `pixi run cargo test -p qrnamap-cli --test io`. Expected: FAIL.
- [ ] **Step 3: Implement `open_input`.** A truncated stream surfaces as `Truncated` when the decoder returns `UnexpectedEof`.
- [ ] **Step 4: Run** the tests. Expected: PASS.
- [ ] **Step 5: Commit** `feat(cli): read plain, gzip and zstd input`.

### Task 5: Transcriptome builder, SAF and splice sites

**Files:**
- Create: `crates/qrnamap-core/src/transcripts.rs`, `crates/qrnamap-core/src/saf.rs`, `crates/qrnamap-core/src/splice.rs`; fixtures in `tests/fixtures/mini/` (`genome.fa`, `genes.gff3`, `README.md`)
- Test: `crates/qrnamap-core/tests/transcripts.rs`, `crates/qrnamap-core/tests/gffread_parity.rs`

**Interfaces:**
- Consumes: `qseq_io::{GeneModel, Gene, Transcript, Strand}`.
- Produces:
  ```rust
  pub struct TxRecord { pub id: String, pub gene: String, pub biotype: String, pub seq: Vec<u8> }
  pub enum BuildError { MissingContig { gene: String, contig: String }, OutOfRange { transcript: String } }
  pub fn build_transcripts(genome: &[(String, Vec<u8>)], model: &GeneModel) -> Result<Vec<TxRecord>, BuildError>
  pub fn write_fasta<W: Write>(records: &[TxRecord], w: W) -> io::Result<()>
  pub fn tx2gene_tsv(records: &[TxRecord]) -> String        // "transcript\tgene\n" rows
  pub fn saf(model: &GeneModel) -> String                   // GeneID Chr Start End Strand; union of exons per gene
  pub fn splice_sites(model: &GeneModel) -> String          // HISAT2 format, 0-based flanking bases, de-duplicated
  ```
- Transcript rule (spec 11c): exons; else CDS merged with explicit UTR features; else CDS only. Minus-strand transcripts are reverse-complemented. rRNA stays in the index; `biotype` carries the `kind`.

- [ ] **Step 1: Make the fixture** by hand: a 3 kb genome with two contigs; genes for one plus-strand two-exon gene, one minus-strand three-exon gene with two isoforms, one rRNA, one CDS-only gene. Document each gene's expected sequence in `README.md`.
- [ ] **Step 2: Write tests** in `transcripts.rs`: `plus_strand_two_exons_concatenate`, `minus_strand_reverse_complemented`, `both_isoforms_emitted`, `cds_only_gene_uses_cds`, `rrna_in_index_with_biotype`, `missing_contig_is_named_error`, `saf_is_exon_union_per_gene` (overlapping exons of two isoforms merge), `tx2gene_has_one_row_per_transcript`.
- [ ] **Step 3: Write `gffread_parity.rs`** (ignored unless `gffread` is on PATH): run `gffread -w` on the fixture and assert every FASTA record equals ours by id and sequence. Add a second case that runs `hisat2_extract_splice_sites.py` (from the hisat2 package) on a GTF made by `gffread -T` and asserts our `splice_sites` output has the same set of lines.
- [ ] **Step 4: Run** `pixi run cargo test -p qrnamap-core`. Expected: FAIL.
- [ ] **Step 5: Implement the four functions** in the three files.
- [ ] **Step 6: Run** `pixi run cargo test -p qrnamap-core -- --include-ignored` and `pixi run wasm-check`. Expected: PASS. A disagreement with `gffread` is reported in the test output, not hidden; decide per case whether our rule or the fixture is wrong and record the decision in `tests/fixtures/mini/README.md`.
- [ ] **Step 7: Commit** `feat(core): transcriptome builder, SAF and splice sites`.

### Task 6: `genome-catalog` crate and generator

**Files:**
- Create: `crates/genome-catalog/src/lib.rs`, `crates/genome-catalog/data/genomes.json` (copied snapshot), `crates/genome-catalog/scripts/make_genome_list.py`, `crates/genome-catalog/scripts/tests/test_make_genome_list.py`
- Test: `crates/genome-catalog/tests/catalog.rs`

**Interfaces:**
- Produces:
  ```rust
  pub struct FileRef { pub name: String, pub url: String, pub bytes: u64 }
  pub struct Genome { pub accession: String, pub source: String, pub label: String, pub organism: String,
      pub strain: Option<String>, pub level: String, pub assembly: String, pub genome_size: u64,
      pub genes: Option<u32>, pub fasta: FileRef, pub gff: Option<FileRef> }
  pub struct Catalog { pub schema_version: u32, pub generated: String, pub genomes: Vec<Genome> }
  pub fn builtin() -> Catalog
  pub fn parse(json: &str) -> Result<Catalog, CatalogError>
  impl Catalog { pub fn find(&self, accession: &str) -> Option<&Genome>;
                 pub fn usable_for_rna(&self) -> Vec<&Genome> }   // gff present
  pub fn size_warning(g: &Genome) -> Option<String>               // Some(...) above 50_000_000 bp
  ```
- The JSON schema is qmutfind's `genomes.json` plus `schema_version`.

- [ ] **Step 1: Copy** qmutfind `bench/make_genome_list.py` and `bench/tests/test_make_genome_list.py` (record the qmutfind commit hash), and a current `site/genomes/genomes.json`; add `"schema_version": 1` to the JSON and to the generator output.
- [ ] **Step 2: Write tests** in `catalog.rs`: `builtin_parses`, `find_known_accession` (`GCF_000002655.1`, Af293), `genome_without_gff_not_usable_for_rna`, `size_warning_above_50mb`, `unknown_schema_version_is_error`. Run the copied pytest with `pixi run pytest crates/genome-catalog/scripts/tests`.
- [ ] **Step 3: Run** `pixi run cargo test -p genome-catalog`. Expected: FAIL.
- [ ] **Step 4: Implement** the types with `serde` and `serde_json`; `builtin()` uses `include_str!`.
- [ ] **Step 5: Run** the Rust tests, the pytest, and `pixi run wasm-check`. Expected: PASS.
- [ ] **Step 6: Commit** `feat(genome-catalog): schema, builtin snapshot, generator from qmutfind@<hash>`.

### Task 7: `genome-fetch` native downloads

**Files:**
- Create: `crates/genome-fetch/src/lib.rs`, `crates/genome-fetch/src/fungidb.rs`, `crates/genome-fetch/src/ncbi.rs`
- Test: `crates/genome-fetch/tests/fetch.rs`

**Interfaces:**
- Consumes: `genome_catalog::{Catalog, Genome}`.
- Produces:
  `pub fn fetch_ncbi(g: &Genome, out_dir: &Path) -> Result<Fetched, FetchError>` (downloads FASTA and GFF3 from the catalog URLs, checks the byte size, writes files as given, already `.gz`)
  `pub fn fetch_fungidb(organism: &str, key: &KeySource, base_url: &str, out_dir: &Path) -> Result<Fetched, FetchError>` (code moved from qmutfind `qmutfind-cli/src/fetch.rs`: key read from `EUPATHDB_KEY` or a key file, sent only to the host of `base_url`, never printed, never over plain http when the site is https)
  `pub struct Fetched { pub fasta: PathBuf, pub gff: Option<PathBuf> }`; `pub enum FetchError { SizeMismatch{..}, NeedsKey{..}, Http{..}, Io(..) }`.
- Ensembl is not built (spec 12).

- [ ] **Step 1: Write tests** against a local server started with `std::net::TcpListener` in the test: `ncbi_download_checks_size` (wrong length → `SizeMismatch`), `ncbi_download_writes_files`, `fungidb_key_sent_only_to_key_host` (redirect to a second host receives no `Authorization` header), `fungidb_401_without_key_names_variable` (error text contains `EUPATHDB_KEY`), `key_never_in_error_text`.
- [ ] **Step 2: Run** `pixi run cargo test -p genome-fetch`. Expected: FAIL.
- [ ] **Step 3: Implement** with `ureq = "2"`. Move the FungiDB code and keep its redirect handling; add the NCBI path as new code.
- [ ] **Step 4: Run** the tests. Expected: PASS.
- [ ] **Step 5: Commit** `feat(genome-fetch): NCBI and FungiDB downloads`.

### Task 8: Sample sheet and design grammar

**Files:**
- Create: `crates/qrnamap-core/src/sheet.rs`, `crates/qrnamap-core/src/design.rs`
- Test: `crates/qrnamap-core/tests/sheet.rs`, `crates/qrnamap-core/tests/design.rs`

**Interfaces:**
- Produces:
  ```rust
  pub enum ColType { Factor, Numeric }
  pub struct Sample { pub name: String, pub fastq1: Vec<String>, pub fastq2: Vec<String>,
      pub covariates: BTreeMap<String, String> }
  pub struct SampleSheet { pub samples: Vec<Sample>, pub columns: Vec<(String, ColType)> }
  pub enum SheetError { DuplicateSample(String), BadName(String), MissingColumn(String),
      MateCountMismatch{sample:String}, BadNumeric{column:String,value:String}, Empty }
  impl SampleSheet { pub fn parse(tsv: &str) -> Result<SampleSheet, Vec<SheetError>> }
  pub struct Design { /* terms over column names */ }
  pub enum DesignError { UnknownColumn(String), BadSyntax(String) }
  impl Design { pub fn parse(formula: &str, columns: &[(String, ColType)]) -> Result<Design, DesignError>;
                pub fn to_r(&self) -> String }
  pub struct Contrast { pub factor: String, pub level: String, pub reference: String }
  ```
- Sheet format (spec 7, my choices for the details): TSV; required columns `sample`, `fastq1`; optional `fastq2`; every other column is a covariate. A header `name:numeric` declares a numeric column; all others are factors. `fastq1` and `fastq2` hold one path or a comma-separated list (lanes). A sample name must match `[A-Za-z][A-Za-z0-9_.-]*`.
- Grammar for the formula: `~` then terms joined by `+`, `:` or `*`; terms are column names. Everything else is `BadSyntax`.

- [ ] **Step 1: Write tests:** sheet — `parses_minimal_sheet`, `lanes_split_on_comma`, `duplicate_sample_is_error`, `name_with_space_rejected`, `name_with_quote_rejected`, `name_with_dollar_rejected`, `mate_lists_must_match_length`, `numeric_header_parses_numbers`, `bad_numeric_value_is_error`, `all_errors_reported_together`. Design — `accepts_condition`, `accepts_batch_plus_condition`, `rejects_unknown_column`, `rejects_function_call` (`~ log(x)`), `rejects_semicolon_and_backtick`, `to_r_quotes_nothing_but_names_and_operators`.
- [ ] **Step 2: Run** `pixi run cargo test -p qrnamap-core --test sheet --test design`. Expected: FAIL.
- [ ] **Step 3: Implement** both files. `to_r` emits only validated column names (already restricted to `[A-Za-z][A-Za-z0-9_.]*`) and the operators `~ + : *`.
- [ ] **Step 4: Run** the tests and `pixi run wasm-check`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(core): sample sheet and design grammar`.

### Task 9: `ToolCall` and native lowering, strandedness map, run config

**Files:**
- Create: `crates/qrnamap-core/src/toolcall.rs`, `crates/qrnamap-core/src/config.rs`
- Test: `crates/qrnamap-core/tests/toolcall.rs`

**Interfaces:**
- Produces:
  ```rust
  pub enum Strandedness { Unstranded, Forward, Reverse }
  pub enum Mode { Quant, Align }
  pub enum Quantifier { Kallisto }                     // Salmon added in Plan 3
  pub struct RunConfig { pub mode: Mode, pub quantifier: Quantifier, pub strandedness: Strandedness,
      pub reads_per_sample: u64,                       // default 20_000_000
      pub fragment_len: Option<(f64, f64)>,            // mean, sd; required for single-end kallisto
      pub max_intron: u32,                             // HISAT2 --max-intronlen, default 5000 (my call; UNVERIFIED for fungi, revisit with data)
      pub threads: u32 }
  pub enum Tool { KallistoIndex, KallistoQuant, Hisat2Build, Hisat2Align, FeatureCounts }
  pub struct ToolCall { pub tool: Tool, pub args: Vec<String>, pub inputs: Vec<PathBuf>, pub outputs: Vec<PathBuf> }
  pub enum ConfigError { FragmentLengthRequired{sample:String} }
  pub fn kallisto_quant(sample:&Sample, index:&Path, out:&Path, cfg:&RunConfig) -> Result<ToolCall, ConfigError>
  pub fn hisat2_align(sample:&Sample, index_prefix:&Path, splice_sites:&Path, sam:&Path, cfg:&RunConfig) -> ToolCall
  pub fn feature_counts(sam:&Path, saf:&Path, out:&Path, paired:bool, cfg:&RunConfig) -> ToolCall
  pub fn argv_native(call:&ToolCall) -> Vec<String>    // program name first
  ```
- Strandedness map (spec 11c): kallisto `--fr-stranded` (Forward) or `--rf-stranded` (Reverse); HISAT2 `--rna-strandness F|R` (single) or `FR|RF` (pairs); featureCounts `-s 1` (Forward), `-s 2` (Reverse), `-s 0` (Unstranded). Always: kallisto `--plaintext -t 1`, HISAT2 `--no-temp-splicesite --known-splicesite-infile`, featureCounts `-p --countReadPairs` when paired and `-F SAF`.

- [ ] **Step 1: Write table tests** (exact argv): `kallisto_paired_unstranded`, `kallisto_paired_reverse_adds_rf_stranded`, `kallisto_single_end_needs_fragment_length` (→ `FragmentLengthRequired`), `kallisto_single_end_uses_l_and_s`, `hisat2_pair_reverse_is_RF`, `hisat2_single_forward_is_F`, `hisat2_has_no_temp_splicesite_and_known_sites`, `featurecounts_paired_has_p_and_countreadpairs`, `featurecounts_reverse_is_s2`, `lanes_are_comma_joined_for_kallisto_and_hisat2`.
- [ ] **Step 2: Run** `pixi run cargo test -p qrnamap-core --test toolcall`. Expected: FAIL.
- [ ] **Step 3: Implement** the types and builders. Lanes: kallisto takes the files in order as arguments; HISAT2 takes `-1 a,b -2 c,d`.
- [ ] **Step 4: Run** the tests and `pixi run wasm-check`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(core): ToolCall, native lowering and strandedness map`.

### Task 10: Count parsers and matrix

**Files:**
- Create: `crates/qrnamap-core/src/counts.rs`
- Test: `crates/qrnamap-core/tests/counts.rs`

**Interfaces:**
- Produces:
  ```rust
  pub fn gene_counts_from_abundance(abundance_tsv:&str, tx2gene:&HashMap<String,String>)
      -> Result<BTreeMap<String,u64>, CountError>        // sum est_counts per gene, round half away from zero
  pub fn parse_featurecounts(tsv:&str) -> Result<BTreeMap<String,u64>, CountError>   // skips the '#' line and header; first column Geneid, last column the count
  pub struct CountMatrix { pub genes: Vec<String>, pub samples: Vec<String>, pub counts: Vec<Vec<u64>> } // counts[gene][sample]
  pub fn assemble(per_sample:&[(String, BTreeMap<String,u64>)], genes:&[String]) -> CountMatrix  // missing gene = 0
  impl CountMatrix { pub fn to_tsv(&self) -> String }
  pub enum CountError { UnknownTranscript(String), BadRow{line:usize} }
  ```

- [ ] **Step 1: Write tests:** `sum_est_counts_per_gene_rounds` (two transcripts 1.4 + 1.2 → 3; 2.5 → 3), `unknown_transcript_is_error`, `featurecounts_skips_comment_and_header`, `featurecounts_reads_last_column`, `assemble_fills_zero`, `matrix_tsv_header_is_gene_then_samples`, `gene_rows_identical_in_both_modes` (genes from `saf` and from `tx2gene` of the Task 5 fixture are the same set).
- [ ] **Step 2: Run** `pixi run cargo test -p qrnamap-core --test counts`. Expected: FAIL.
- [ ] **Step 3: Implement** the functions.
- [ ] **Step 4: Run** the tests. Expected: PASS.
- [ ] **Step 5: Commit** `feat(core): count parsers and gene-by-sample matrix`.

### Task 11: Native pipelines and strandedness check (`transcripts`, `quant`, `align`, `matrix`, `strand-check`)

**Files:**
- Create: `crates/qrnamap-cli/src/exec.rs`, `crates/qrnamap-cli/src/pipeline.rs`, `crates/qrnamap-cli/src/main.rs` (subcommands), `crates/qrnamap-core/src/strand.rs`; a read simulator in `crates/qrnamap-cli/tests/common/sim.rs`
- Test: `crates/qrnamap-cli/tests/pipeline_quant.rs`, `crates/qrnamap-cli/tests/pipeline_align.rs`, `crates/qrnamap-core/tests/strand.rs`

**Interfaces:**
- Consumes: Tasks 2 to 10.
- Produces: `pub fn run_call(call:&ToolCall, scratch:&Path) -> Result<(), ExecError>` (runs `argv_native`, checks exit code and that every output exists, stderr in the error); CLI subcommands exactly as in spec 8 except `fetch`/`de`/`run` (Tasks 12 and 13): `qrnamap transcripts`, `quant`, `align`, `matrix`, `strand-check`. Output tables are written `.zst` unless the path ends `.gz` or is plain. Temporary files use `${SCRATCH:?}` when set, otherwise the `--tmpdir` flag.
  `pub fn assigned_fraction(strand_counts: (u64,u64,u64)) -> StrandReport` in `strand.rs`: `(sense, antisense, ambiguous)` → fractions and a `suggested: Strandedness`; warn when the winning fraction is below 0.8 (my threshold, UNMEASURED).
- The simulator writes error-free paired FASTQ from a `Vec<TxRecord>` and a vector of read counts per transcript, so truth is known.

- [ ] **Step 1: Write tests:** `strand.rs` — `forward_library_suggests_forward`, `unstranded_library_warns_low_fraction`. Pipeline (skipped when the tool is missing, run in CI under pixi) — `quant_recovers_simulated_counts_within_10_percent` (10,000 pairs over the mini fixture transcripts; gene counts within 10% of truth, tolerance is my choice, tighten after measuring), `align_counts_match_truth_within_5_percent` (HISAT2 + featureCounts on reads simulated from the genome including a spliced gene), `both_modes_same_gene_rows`, `quant_with_lanes_equals_concatenated`, `native_failure_reports_stderr` (a bad flag → `ExecError` containing tool output).
- [ ] **Step 2: Run** `pixi run cargo test -p qrnamap-cli -p qrnamap-core --test strand`. Expected: FAIL.
- [ ] **Step 3: Implement** `exec.rs`, `pipeline.rs` and the subcommands. Each sample is an independent function `quant_sample` / `align_sample` so Task 13 can wrap them in SLURM tasks. `hisat2-build` reads a decompressed FASTA in `$SCRATCH` (reading gzip directly is UNVERIFIED). `strand-check` maps the first 100,000 pairs (parameter) with the strand-aware featureCounts output and prints `StrandReport`; the same report is printed after `align` and `quant` runs as "assigned fraction" per sample (for `quant`, from `run_info.json` `p_pseudoaligned`).
- [ ] **Step 4: Run** `pixi run cargo test --workspace`. Expected: PASS. Record the measured recovery errors in `README.md` under "Measured".
- [ ] **Step 5: Commit** `feat(cli): quant, align, matrix and strand-check`.

### Task 12: DESeq2 step (`de`)

**Files:**
- Create: `crates/qrnamap-core/src/de.rs`, `crates/qrnamap-cli/src/de.rs`
- Test: `crates/qrnamap-core/tests/de_script.rs`, `crates/qrnamap-cli/tests/de_run.rs`

**Interfaces:**
- Consumes: `SampleSheet`, `Design`, `Contrast`, `CountMatrix`.
- Produces:
  ```rust
  pub fn check_design(sheet:&SampleSheet, contrast:&Contrast) -> Result<(), DeError>  // >=2 samples per tested group, no all-zero sample (given the matrix)
  pub fn deseq2_script(sheet:&SampleSheet, design:&Design, contrast:&Contrast, counts_path:&str, out_prefix:&str) -> String
  pub enum DeError { TooFewReplicates{group:String,n:usize}, AllZeroSample(String), NotFullRank(String) }
  ```
  CLI: `qrnamap de --counts counts.tsv --samples sheet.tsv --design '~ condition' --contrast condition,treated,control [--shrink normal|apeglm|ashr]`. The script writes `<prefix>.results.csv`, `<prefix>.normalized.csv`, `<prefix>.vst.csv`, and prints `packageVersion("DESeq2")` and `R.version.string` into `<prefix>.versions.txt`. Default: Wald test, no shrinkage. Uses `varianceStabilizingTransformation` (not `vst()`, spec 11c). Sample and level names are written as R string literals with `"` and `\` escaped (names are already restricted by Task 8; this is a second guard). A rank-deficient design shows R's error text in `DeError::NotFullRank`.

- [ ] **Step 1: Write tests:** `script_default_has_no_shrink`, `script_uses_varianceStabilizingTransformation`, `script_sets_reference_level_explicitly`, `numeric_column_not_converted_to_factor`, `one_replicate_group_refused`, `all_zero_sample_refused`, `script_contains_no_unquoted_user_text` (feed names with `"` and `\`, which Task 8 normally rejects, built by hand; assert escaping). Run test (`de_run.rs`, skipped without Rscript): a 200-gene synthetic matrix with 3 vs 3 samples and 20 genes truly changed 4-fold → at least 15 of the 20 have `padj < 0.05` and no more than 5 of the other 180 do (thresholds are my choice for a synthetic case). Also `confounded_batch_is_not_full_rank` (batch equals condition → `NotFullRank`).
- [ ] **Step 2: Run** `pixi run cargo test -p qrnamap-core --test de_script`. Expected: FAIL.
- [ ] **Step 3: Implement** `check_design`, `deseq2_script`, and the CLI wrapper that runs `Rscript` and maps R failures to `DeError`.
- [ ] **Step 4: Run** `pixi run cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat: DESeq2 step with design checks`.

### Task 13: `fetch`, `run`, SLURM array, `reads_per_sample` check

**Files:**
- Create: `crates/qrnamap-cli/src/fetch_cmd.rs`, `crates/qrnamap-cli/src/run.rs`, `crates/qrnamap-cli/src/slurm.rs`
- Test: `crates/qrnamap-cli/tests/run.rs`, `crates/qrnamap-cli/tests/slurm.rs`

**Interfaces:**
- Produces: `qrnamap fetch --source ncbi|fungidb --id <accession-or-organism> --out DIR [--eupathdb-key-file F]`; `qrnamap run --samples S --reference R --gff G --mode quant|align [flags of Tasks 11 and 12]`;
  `pub fn count_pairs(path:&Path) -> Result<u64, InputError>` (uses `open_input` and `FastqReader`);
  `pub fn check_reads_per_sample(sheet:&SampleSheet, limit:u64) -> Vec<Warning>`;
  `pub fn slurm_array_script(plan:&Plan) -> String` where `Plan { samples: Vec<String>, per_task: usize, cpus: u32, mem_gb: u32, time: String, run_dir: PathBuf }`.
- `run` order: build transcripts and index once (index on shared storage), then per-sample steps, then `matrix`, then `de` if a contrast is given. `--slurm` writes the array script and an `afterok` job for `matrix` and `de`, and does not run tools. The script reads `$SCRATCH` with `${SCRATCH:?}`, takes paths from arguments or the environment, and never uses `BASH_SOURCE[0]`. Samples per task: the CLI flag `--samples-per-task`; the default is 1 until a measured per-sample time exists (global rule: size tasks to about 1 to 1.5 hours from real timing, not a guess).

- [ ] **Step 1: Write tests:** `count_pairs_plain_and_gzip`, `exceeds_reads_per_sample_gives_warning_not_error`, `slurm_script_uses_scratch_with_required_check` (contains `${SCRATCH:?}`), `slurm_script_has_no_bash_source` (does not contain `BASH_SOURCE`), `slurm_array_range_matches_task_count` (7 samples, 3 per task → `--array=0-2`), `slurm_index_job_precedes_array` (array depends on the index job with `afterok`), `fetch_ncbi_writes_files_to_out_dir` (local server, reusing Task 7's test helper), `run_end_to_end_quant_on_mini_fixture` (skipped without kallisto; produces `counts.tsv` and a `results.csv`).
- [ ] **Step 2: Run** `pixi run cargo test -p qrnamap-cli --test run --test slurm`. Expected: FAIL.
- [ ] **Step 3: Implement** the three modules and wire them into `main.rs`.
- [ ] **Step 4: Run** `pixi run cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(cli): fetch, run and SLURM array`.

### Task 14: Real-data validation on the cluster (CLI side of spec step 7)

**Files:**
- Create: `validation/README.md`, `validation/run_yeast.sh`, `validation/compare_gffread.py`
- Output (not committed): `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/`

**Interfaces:**
- Consumes: the full CLI. Produces: a "Measured" section in `README.md` with numbers from these runs only.

- [ ] **Step 1: Choose and record the data.** Pick a small public yeast RNA-seq set with at least 2 conditions and 3 replicates each from the ENA portal API (`fields=read_count,base_count,library_layout`); write accessions, read counts and layout in `validation/README.md`. Download with `prefetch` + `parallel-fastq-dump` (global rule; do not use ENA FTP) and check both mates equal the ENA `read_count`. Fetch the *S. cerevisiae* S288C genome (`GCF_000146045.2`) with `qrnamap fetch`.
- [ ] **Step 2: Transcript parity.** Run `qrnamap transcripts` and `gffread -w` on S288C; `compare_gffread.py` (reads `.gz`/`.zst`) lists records that differ. Expected: zero differences, or each difference explained in `validation/README.md`.
- [ ] **Step 3: Both modes.** Run `qrnamap run --mode quant` and `--mode align` with `--slurm` (3 samples per task is a guess; record real per-sample time and set `--samples-per-task` from it). Check that both gene matrices have identical row names. Record per-sample mapping rate, assigned fraction, and wall time.
- [ ] **Step 4: Agreement.** Compare per-gene counts between modes (Spearman correlation) and DE results between modes (overlap of genes with `padj < 0.05`). Report the numbers as measured; do not describe them as good or bad without a stated criterion.
- [ ] **Step 5: Strandedness.** Run `strand-check` on every sample; confirm its suggestion matches the library kit stated in the ENA record, or record that the kit is unknown.
- [ ] **Step 6: Commit** `docs: CLI validation on yeast RNA-seq (measured)`.

---

## Self-Review (run 2026-11-19)

- **Spec coverage:** goals and pipeline (spec 1, 6): Tasks 5, 10, 11, 12. Sample sheet (7): Task 8. CLI (8): Tasks 11 to 13. Compressed reads (11f): Task 4. Crates (5): Tasks 2 to 7. Catalog (5, 11a): Task 6. Corrections 11c (all isoforms, SAF, strandedness table, design grammar, DE version note, `--known-splicesite-infile`, gffread in pixi, parity with `-t 1`): Tasks 3, 5, 8, 9, 11, 12. Reads parameter and licence (11a): Task 9 and `LICENSE`. Validation CLI side (10 step 7): Task 14. **Not covered here on purpose:** spike 0 and everything browser (Plans 2 and 3); Salmon; chunked streaming; `lengthScaledTPM`; DESeq2 1.52.0 swap; Ensembl.
- **Gaps found and fixed:** the plan first used `FastqChunker` for chunking; chunking is browser work, so only the reader is used here (Task 2 keeps the chunker for Plan 2).
- **Type consistency:** `GeneModel`/`Gene`/`Transcript` (Task 3) are the inputs of Task 5; `Sample`, `SampleSheet`, `Design`, `Contrast` (Task 8) are used in Tasks 9, 11, 12; `RunConfig`, `ToolCall` (Task 9) in Tasks 11 and 13; `CountMatrix` (Task 10) in Tasks 11 and 12.
- **Values I chose and the spec did not fix (all open to change):** `max_intron` 5000; strand-fraction warning threshold 0.8; simulation tolerances (10%, 5%); synthetic DE thresholds; `.zst` default for output tables; first-100,000-pair strand check; sheet header `name:numeric`.
