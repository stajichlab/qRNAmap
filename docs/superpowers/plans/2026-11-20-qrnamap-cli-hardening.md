# qRNAmap CLI Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the open items of the CLI-and-core plan that the user chose to do first: write gene-class labels and per-sample mapping summaries as sidecar files, lower and document the SLURM memory default from the measured memory behaviour, let `run --slurm` work without a site-set `$SCRATCH`, and fold the AMB experiment results into the docs.

**Architecture:** The count table `counts.tsv.zst` is unchanged (DESeq2 reads it). Two new sidecar files sit next to it: `genes.tsv.zst` (one row per gene: class, biotype, mitochondrial flag) and `mapping_summary.tsv.zst` (one row per sample). Gene classes are computed in `qrnamap-core` (pure); files are written by `qrnamap-cli` through the existing atomic `Output` writer. SLURM script changes stay inside `slurm.rs`.

**Tech Stack:** Rust (the existing workspace), kallisto 0.52.0, HISAT2 2.2.2, featureCounts, R/DESeq2 (pixi environment), SLURM scripts.

**Spec:** `docs/superpowers/specs/2026-11-19-qrnamap-design.md` (sections 6, 8, 11a, 11c). Evidence for the memory and scratch tasks: `validation/amb/README.md` (experiments E1 to E3 and defect D3) and `validation/README.md`. User decisions recorded 2026-11-20: sidecar files with the matrix unchanged (option A); default `--mem-gb` 8 (option A, default taken because no answer had arrived; the user may change it).

## Global Constraints

- Licence GPL-3.0-only. pixi manifests use `[workspace]`, never `[project]`.
- `qseq-io`, `genome-catalog`, `qrnamap-core` have no file or network I/O and must build with `cargo build -p qrnamap-core -p qseq-io -p genome-catalog --target wasm32-unknown-unknown`.
- Large text output is compressed by the output path's suffix (`.zst`, `.gz`, plain); new CLI tables default to `.zst` names. Writers use the existing atomic `.partial` + rename path.
- Ruling R24: every step removes the outputs a previous run may have left BEFORE it runs any tool, so a failed rerun cannot leave stale files that look complete.
- Ruling R18: paths and sample names reach tools only as argv; paths baked into shell text go through `sh_quote`; no sample names or FASTQ paths in any generated script.
- Global lab rules: `${SCRATCH:?}` for node-local scratch where `--tmpdir` is not given; never hardcode `/scratch/...`; never `BASH_SOURCE`; `--time`, `--mem-gb`, `--cpus` defaults are labelled UNMEASURED placeholders unless a measurement is cited.
- Ruling R5: tests needing an external tool skip with a message when the tool is missing, and FAIL under `QRNAMAP_REQUIRE_TOOLS=1`. Ruling R11: run tools as `PATH=<worktree>/.pixi/envs/default/bin:$PATH ...`; never `pixi run` (it hangs on the UCR login node).
- Gene labels containing a tab, newline or CR are refused before any table is written (existing rule; it also applies to `genes.tsv.zst`).
- Commit messages end with `Co-Authored-By: <model> <noreply@anthropic.com>` and the line `Claude-Session: https://claude.ai/code/session_01CoSBcPmuXVDpBrrPNsnYtc`.
- Never report a number that was not measured. Do not describe a number as good or bad without a stated criterion.

## Review Focus

1. A GFF3 without `gene_biotype` and with unusual transcript types (`transcript`, `primary_transcript`, `lnc_RNA`, `V_gene_segment`): the gene still gets a row, class `other`, no panic. (Task 2)
2. A rerun that fails leaves an old `genes.tsv.zst` or `<sample>.summary.tsv.zst`: both are removed first, and `matrix` refuses a missing or foreign summary with a named error. (Tasks 3, 4)
3. `--mito-contig` names a contig that has no genes, or the GFF has a `genome=mitochondrion` region for a contig with no genes: no error, no empty row. (Tasks 1, 2)
4. A `--tmpdir` or `--out` with a space, `$`, quote or backtick in `run --slurm`: scripts parse (`bash -n`) and the path arrives as one literal argument. (Task 5)
5. A job submitted from a bare environment (no `$SCRATCH`): the job stops at its first line with a message that names the cause and both remedies. (Task 5)

## File Structure

```
crates/qseq-io/src/gff.rs                  Gene.biotype, GeneModel.mito_contigs (modify)
crates/qrnamap-core/src/genes.rs           GeneClass, GeneRow, gene_table, gene_table_tsv (create)
crates/qrnamap-core/src/lib.rs             re-export (modify)
crates/qrnamap-core/tests/genes.rs         core tests (create)
crates/qrnamap-cli/src/pipeline.rs         write_genes_table, write_sample_summary, build_mapping_summary (modify)
crates/qrnamap-cli/src/run.rs, main.rs     wire flags and outputs, stale-output removal (modify)
crates/qrnamap-cli/src/slurm.rs            default mem, --tmpdir, SCRATCH message (modify)
crates/qrnamap-cli/tests/*.rs              new CLI tests (create/modify)
README.md, validation/README.md, validation/amb/*   docs (modify)
```

---

### Task 1: The gene model keeps `gene_biotype` and mitochondrial regions

**Files:**
- Modify: `crates/qseq-io/src/gff.rs`
- Test: `crates/qseq-io/tests/gff.rs`

**Interfaces:**
- Consumes: the existing `parse_gff3`, `Gene`, `GeneModel`, `Feat.attr(key)`.
- Produces: `Gene.biotype: Option<String>` (the value of the `gene_biotype` attribute of the gene-type feature, percent-decoded, `None` when absent or empty); `GeneModel.mito_contigs: std::collections::BTreeSet<String>` (the seqid of every `region` feature whose `genome` attribute equals `mitochondrion` or `mitochondrial`, case-sensitive). `GeneModel` derives `Default`.

- [ ] **Step 1: Write the failing tests** in `tests/gff.rs`: `gene_biotype_attribute_is_kept` (a `gene` with `gene_biotype=rRNA` gives `genes[0].biotype == Some("rRNA".to_string())`); `missing_gene_biotype_is_none`; `percent_encoded_biotype_is_decoded`; `mito_region_contig_is_collected` (a `region` line `NC_001224.1 . region 1 85779 . + . ID=NC_001224.1:1..85779;genome=mitochondrion;Name=MT` puts `NC_001224.1` in `mito_contigs`); `chromosome_region_is_not_mito` (`genome=chromosome`); `mito_region_without_genes_is_not_an_error`.
- [ ] **Step 2: Run** `PATH=<worktree>/.pixi/envs/default/bin:$PATH cargo test -p qseq-io --test gff`. Expected: FAIL (fields missing).
- [ ] **Step 3: Implement** the two fields in `gff.rs`; capture `gene_biotype` where a gene is first built; collect region seqids in the pass-1 feature loop (regions are not children and not gene-type, so add one branch). Update every place that builds `Gene` or `GeneModel` by hand in the workspace (the compiler lists them; use `biotype: None` and `Default::default()`).
- [ ] **Step 4: Run** `cargo test --workspace` and the wasm32 build of the three pure crates. Expected: PASS, no existing assertion changed.
- [ ] **Step 5: Commit** `feat(qseq-io): keep gene_biotype and mitochondrial regions in the gene model`.

### Task 2: The gene table in `qrnamap-core`

**Files:**
- Create: `crates/qrnamap-core/src/genes.rs`, `crates/qrnamap-core/tests/genes.rs`
- Modify: `crates/qrnamap-core/src/lib.rs`

**Interfaces:**
- Consumes: `qseq_io::{GeneModel, Gene, Strand}` from Task 1.
- Produces:
  ```rust
  pub enum GeneClass { ProteinCoding, Rrna, Trna, OtherNcrna, Pseudogene, Other }
  impl GeneClass { pub fn as_str(&self) -> &'static str }   // "protein_coding" "rRNA" "tRNA" "other_ncRNA" "pseudogene" "other"
  pub struct GeneRow { pub gene: String, pub contig: String, pub strand: char, pub biotype: String,
                       pub class: GeneClass, pub mito: bool, pub n_transcripts: usize }
  pub fn gene_table(model: &GeneModel, extra_mito_contigs: &[String]) -> Vec<GeneRow>   // model order
  pub fn gene_table_tsv(rows: &[GeneRow]) -> String   // header `gene\tcontig\tstrand\tbiotype\tclass\tmito\tn_transcripts`, `mito` is `yes`/`no`
  ```
- Rule (exact): `biotype` = `gene.biotype`, else `"mRNA"` when any transcript has kind `mRNA`, else the first transcript's kind, else `"unknown"`. Class from the lower-cased biotype: `protein_coding` or `mrna` -> ProteinCoding; `rrna` -> Rrna; `trna` -> Trna; contains `pseudogene` or `pseudogenic` -> Pseudogene; ends with `rna` (covers snoRNA, snRNA, ncRNA, lnc_RNA, antisense_RNA, misc_RNA, SRP_RNA, RNase_P_RNA, RNase_MRP_RNA, telomerase_RNA) -> OtherNcrna; anything else -> Other. `mito` is true when the gene's contig is in `model.mito_contigs` or in `extra_mito_contigs`.

- [ ] **Step 1: Write the failing tests** in `tests/genes.rs` with literal expectations: `classes_for_the_yeast_biotypes` (a table of biotype -> class for `protein_coding`, `rRNA`, `tRNA`, `snoRNA`, `snRNA`, `ncRNA`, `antisense_RNA`, `misc_RNA`, `RNase_P_RNA`, `SRP_RNA`, `telomerase_RNA`, `pseudogene`, `lnc_RNA`, `V_gene_segment`=other); `biotype_falls_back_to_transcript_kind`; `mrna_wins_over_other_transcript_kinds`; `unknown_when_no_biotype_and_no_transcripts`; `mito_from_region_and_from_extra_contig`; `extra_mito_contig_without_genes_is_fine`; `tsv_header_and_rows_exact` (full literal string for a 3-gene model); `rows_follow_model_order`.
- [ ] **Step 2: Run** `cargo test -p qrnamap-core --test genes`. Expected: FAIL.
- [ ] **Step 3: Implement** `genes.rs` and re-export it from `lib.rs`.
- [ ] **Step 4: Run** the core tests and the wasm32 build. Expected: PASS.
- [ ] **Step 5: Commit** `feat(core): gene table with class, biotype and mitochondrial flag`.

### Task 3: The CLI writes `genes.tsv.zst`

**Files:**
- Modify: `crates/qrnamap-cli/src/pipeline.rs`, `crates/qrnamap-cli/src/run.rs`, `crates/qrnamap-cli/src/main.rs`
- Test: `crates/qrnamap-cli/tests/genes_table.rs` (create)

**Interfaces:**
- Consumes: `qrnamap_core::{gene_table, gene_table_tsv}`, the existing `Output` writer, `prepare_quant`, `prepare_align`, `remove_old_run_outputs` (names as in the code).
- Produces: `pub fn write_genes_table(path: &Path, model: &GeneModel, mito_contigs: &[String]) -> Result<(), PipelineError>` (atomic, suffix decides compression); flag `--mito-contig NAME` (repeatable) on `transcripts`, `quant`, `align` and `run`; `transcripts --genes-table PATH` (optional); `<out>/genes.tsv.zst` written by `prepare_quant`, `prepare_align` and `run --step prepare`; `genes.tsv.zst` added to the files removed at the start of prepare and of plain `quant`/`align` (R24).

- [ ] **Step 1: Write the failing tests** (mini fixture, tool-gated where a tool runs): `transcripts_writes_genes_table` (the mini fixture's rRNA gene has class `rRNA`, the CDS-only gene `protein_coding` or per its kind as the Task 2 rule says; assert exact rows read back with `zstd`); `mito_contig_flag_marks_genes`; `plain_quant_writes_genes_table_in_out` and `run_prepare_writes_genes_table`; `old_genes_table_is_removed_before_a_failing_run` (leave a stale file, run with a missing FASTQ, assert the file is gone); `genes_table_label_with_tab_is_refused`; `counts_matrix_is_unchanged` (the first line of `counts.tsv.zst` is still `gene\t<sample>...`).
- [ ] **Step 2: Run** `cargo test -p qrnamap-cli --test genes_table`. Expected: FAIL.
- [ ] **Step 3: Implement** `write_genes_table`, the flag plumbing, and the stale-output removal.
- [ ] **Step 4: Run** `PATH=<worktree>/.pixi/envs/default/bin:$PATH QRNAMAP_REQUIRE_TOOLS=1 cargo test --workspace`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(cli): write genes.tsv.zst next to the count matrix`.

### Task 4: Per-sample summary files and `mapping_summary.tsv.zst`

**Files:**
- Modify: `crates/qrnamap-cli/src/pipeline.rs`, `crates/qrnamap-cli/src/main.rs`
- Test: `crates/qrnamap-cli/tests/mapping_summary.rs` (create)

**Interfaces:**
- Consumes: `SampleResult.summary: SampleSummary { reads_processed, fraction_mapped, fraction_assigned }`, `QuantContext`/`AlignContext` config (`Strandedness`), `remove_old_outputs(out_dir, sample)`, the `strand_check.tsv` file written by `strand-check` (read its format in the code).
- Produces: `<out>/<sample>.summary.tsv.zst` written atomically by `quant_sample` and `align_sample` as the last step (header `sample\tmode\treads_processed\tfraction_mapped\tfraction_assigned\tstrandedness\tqrnamap_version`, one row); removed by `remove_old_outputs` before any tool runs; `pub fn build_mapping_summary(out_dir: &Path, samples: &[String], strand_check: Option<&Path>) -> Result<String, PipelineError>` (rows in sheet order, extra last column `strand_check_suggestion`, empty when no `strand_check.tsv` row exists for the sample); `qrnamap matrix` writes `<out>/mapping_summary.tsv.zst`; a missing per-sample summary or one whose `sample` cell differs from the expected name is a named error that tells the user to rerun that sample.

- [ ] **Step 1: Write the failing tests:** `quant_sample_writes_summary_equal_to_printed_numbers` and `align_sample_writes_summary_equal_to_printed_numbers` (tool-gated, mini or large fixture; assert the values equal `SampleResult.summary`); `summary_removed_before_a_failing_rerun`; `matrix_writes_mapping_summary_in_sheet_order`; `matrix_refuses_missing_summary_naming_the_sample`; `matrix_refuses_summary_for_another_sample`; `strand_check_suggestion_is_joined_when_present` and `..._empty_when_absent`; `summary_header_exact`.
- [ ] **Step 2: Run** `cargo test -p qrnamap-cli --test mapping_summary`. Expected: FAIL.
- [ ] **Step 3: Implement** the writer, the cleanup and `build_mapping_summary`; wire it into `matrix`.
- [ ] **Step 4: Run** the full workspace tests with tools required. Expected: PASS; no existing assertion changed.
- [ ] **Step 5: Commit** `feat(cli): per-sample summaries and mapping_summary.tsv.zst`.

### Task 5: SLURM memory default, `--tmpdir` for `run --slurm`, clearer scratch message

**Files:**
- Modify: `crates/qrnamap-cli/src/slurm.rs`, `crates/qrnamap-cli/src/run.rs`, `crates/qrnamap-cli/src/main.rs`
- Test: `crates/qrnamap-cli/tests/slurm.rs`

**Interfaces:**
- Consumes: `Plan`, `sh_quote`, the generated-script writers, `DEFAULT_MEM_GB`.
- Produces: `DEFAULT_MEM_GB = 8` (comment: placeholder; measured on AMB data: HISAT2 41,208 KB and featureCounts 36,144 KB resident, the 13 GB in `sacct` was page cache of scratch files, a `--mem 4G` job finished with identical counts; not measured for genomes near 50 Mb or 20 million pairs); `Plan.tmpdir: Option<PathBuf>`; with `--tmpdir DIR` on `run --slurm` the scripts pass `--tmpdir <sh_quote(DIR)>` and contain no `${SCRATCH:?}` check; without it behaviour is as today but the failing message reads: `SCRATCH is not set: this job needs node-local scratch. Submit from a login shell where $SCRATCH is set, or regenerate the scripts with --tmpdir DIR.`; the help text for `--mem-gb` states the default is an UNMEASURED placeholder and why it is small.

- [ ] **Step 1: Write the failing tests:** `default_mem_is_8_gb_in_scripts` (`#SBATCH --mem=8G`); `tmpdir_option_replaces_scratch_check`; `without_tmpdir_scratch_check_remains_with_new_message` (exact message text); `hostile_tmpdir_round_trips` (space, `$HOME`, `'`, backtick: `bash -n` passes and the stub-`sbatch`/stub-`qrnamap` run shows the directory as one literal argument); `tmpdir_with_newline_is_refused`; `mem_gb_zero_is_still_an_error`; `slurm_scripts_never_contain_sample_names_or_fastq_paths` still holds; update the existing tests that asserted 16 G.
- [ ] **Step 2: Run** `cargo test -p qrnamap-cli --test slurm`. Expected: FAIL.
- [ ] **Step 3: Implement** the default, `Plan.tmpdir`, the flag on `run`, and the message.
- [ ] **Step 4: Run** the full workspace tests with tools required and clippy (`cargo clippy --workspace --all-targets`: 0 warnings). Expected: PASS.
- [ ] **Step 5: Commit** `feat(slurm): 8 GB default memory, --tmpdir for run --slurm, clearer SCRATCH message`.

### Task 6: Docs and the AMB results

**Files:**
- Modify: `README.md`, `validation/README.md`, `validation/amb/README.md`, `validation/amb/results/e3_memory.tsv`, `docs/superpowers/execution/2026-11-19-qrnamap-cli-core-ledger.md`

**Interfaces:** none (documentation). Every number must come from `validation/README.md` or `validation/amb/README.md` or from the code.

- [ ] **Step 1: README.md:** document `genes.tsv.zst` (columns, class rule, the GFF `gene_biotype` and `genome=mitochondrion` sources, `--mito-contig`), `mapping_summary.tsv.zst` (columns), `run --slurm --tmpdir`, the 8 GB default and the memory explanation (page cache; the three measured values), the D3 note (submit from a login shell or use `--tmpdir`), the `est_counts` rounding note (gene count = sum of `est_counts` per gene rounded once, half away from zero), and a pointer to `validation/amb/README.md` with its headline numbers copied verbatim. Update "Not in this version".
- [ ] **Step 2: validation/README.md:** add an index line to `validation/amb/README.md` and state that the D1 fix (f2bfb54) was run on real jobs in the AMB experiments (E1, E2) with the outcome as written there.
- [ ] **Step 3: validation/amb/README.md:** soften the overclaims the review found: D3 job 3 ran on one node (r47), not "this cluster"; replace "none of my 25 jobs is queued" by "none of the AMB jobs"; list the two failed probe jobs 32309811 and 32309812 (a reviewer's earlier D3 attempts, output path unwritable); fix the stale header row of `results/e3_memory.tsv` ("all three E3 jobs r46" -> per-row node column only).
- [ ] **Step 4: Ledger:** append the rulings of this plan (R30: sidecar files, user option A; R31: `--mem-gb` default 8; R32: D3 handled by `--tmpdir` plus a clearer message, no site-specific `module load`) with their costs if wrong.
- [ ] **Step 5: Verify** `cargo test --workspace` still passes (docs only, quick) and commit `docs: sidecar files, memory and scratch notes, AMB results index`.

---

## Self-Review (2026-11-20)

- **Spec coverage:** spec 11c "flag these genes in the matrix" -> Tasks 1 to 3 (as a sidecar, per the user's choice); spec 2 "mapping summaries" -> Task 4; spec 6 "stated in the output notes" (rounding) -> Task 6 (README); spec 11a resource requests -> only the memory default (Task 5); time from `reads_per_sample` is NOT in this plan (open, user). Open items not here on purpose: apeglm/ashr packages, the browser plan, Salmon, the comparison with the authors' published GSE80357 tables (an experiment; asked of the user separately), push of `main`.
- **Type consistency:** `GeneModel.mito_contigs` and `Gene.biotype` (Task 1) are consumed by `gene_table` (Task 2) and `write_genes_table` (Task 3); `SampleSummary` fields (existing) feed `write_sample_summary` and `build_mapping_summary` (Task 4); `Plan.tmpdir` (Task 5) is the only new `Plan` field.
- **Values I chose where the spec is silent (open to change):** class names and the `ends_with("rna")` rule; `mito` as `yes`/`no`; `<sample>.summary.tsv.zst` file name; `strand_check_suggestion` last column; the 8 GB default; no `module load` in scripts.
