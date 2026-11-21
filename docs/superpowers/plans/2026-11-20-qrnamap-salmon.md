# Salmon quantifier support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `--quantifier kallisto|salmon` to `qrnamap quant` and `run`, so Salmon output feeds the same count, gene, summary and DE files as kallisto.

**Architecture:** Branch on `Quantifier` inside the existing code. Pure Salmon pieces (argv builders, `quant.sf` gene sum) go in `qrnamap-core`. The CLI crate picks the tool in `prepare_quant`, `quant_sample`, the marker and the SLURM plan. Kallisto code paths do not change.

**Tech Stack:** Rust workspace, Salmon 2.8.0 (bioconda, via pixi), existing test fixtures in `crates/qrnamap-cli/tests/large.rs`.

**Spec:** `docs/superpowers/specs/2026-11-20-salmon-support-design.md` (read it first, including its last section "Decisions after the builder check"). Transcripts-only index. No decoys, no WASM, no align-mode Salmon.

## Global Constraints

- Default quantifier is `kallisto`. Kallisto behavior and output (counts, gene counts) must not change. Only the summary file format changes (new column, Task 2).
- Run tools with `PATH=<worktree>/.pixi/envs/default/bin:$PATH`. Never `pixi run` (ruling R11: it hangs on the login node). Use `CARGO_TARGET_DIR=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-harden`.
- Pure crates (`qseq-io`, `genome-catalog`, `qrnamap-core`) must still build for `wasm32-unknown-unknown`. No JSON dependency in core.
- R18: argv only, no user text in a shell. R24: remove stale outputs first (check inputs, remove, build). R27: marker mismatch refuses `--task-id` and `--step final`. R5: tool tests skip without the tool and fail under `QRNAMAP_REQUIRE_TOOLS=1`.
- Salmon library types: paired `IU`/`ISF`/`ISR`, single-end `U`/`SF`/`SR` for unstranded/forward/reverse. Never `-l A`.
- Salmon paired lanes: all R1 files after `-1`, all R2 files after `-2`. Single-end: all files after `-r`. Kallisto uses adjacent pairs; do not touch that.
- Gene counts: sum `NumReads` per gene through tx2gene, round once, half away from zero. Negative, NaN or infinite `NumReads` is `BadRow`.
- Commit trailers on every commit: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01CoSBcPmuXVDpBrrPNsnYtc`.
- Writing style in docs and comments: Simplified Technical English, short sentences, no inflated claims. Report measured numbers only.

## Review Focus

- Salmon run with `--task-id` on a directory prepared with kallisto: must refuse with a message naming the quantifier (Task 4).
- Switching tools in the same output directory: no stale `abundance.tsv`, `quant.sf` or old summary survives (Task 3).
- Zero mapped reads with Salmon: named `NoReadsPseudoaligned` error, not a silent zero sample (Task 3).
- Two lanes of paired reads give the same counts as the concatenated reads (Task 5).
- A 7-column summary file left by an old version: clear error, no panic (Task 2).
- `--quantifier salmon` with `--mode align`: error, not silently ignored (Task 4).

## File Structure

- Modify `pixi.toml`.
- Modify `crates/qrnamap-core/src/config.rs`, `toolcall.rs`, `counts.rs` and their tests in `crates/qrnamap-core/tests/{toolcall,counts}.rs`.
- Modify `crates/qrnamap-cli/src/pipeline.rs` (summary, `prepare_quant`, `quant_sample`, errors, removal), `run.rs` (marker, `prepared_quant`, `required_tools`, `slurm_plan`), `main.rs` (flag).
- Modify tests: `mapping_summary.rs`, `cli.rs`, `run.rs`, `slurm.rs`, `pipeline_quant.rs`, `large.rs`.
- Modify docs: `README.md`, `validation/README.md`.

---

### Task 1: Environment, config and pure core pieces

**Files:**
- Modify: `pixi.toml`, `crates/qrnamap-core/src/{config,toolcall,counts}.rs`
- Test: `crates/qrnamap-core/tests/{toolcall,counts}.rs`

**Interfaces:**
- Produces (used by Tasks 2 to 5):
  - `Quantifier::Salmon`, and `impl Quantifier { pub fn name(self) -> &'static str }` returning `"kallisto"` or `"salmon"`.
  - `Tool::SalmonIndex`, `Tool::SalmonQuant` (program name `salmon`).
  - `pub fn salmon_index(transcripts_fasta: &str, index_dir: &str, cfg: &RunConfig) -> Result<ToolCall, ConfigError>`. Args: `index -t <fa> -i <dir> -p <threads> --keepDuplicates`. `outputs` lists the directory. (`--keepDuplicates` keeps sequence-identical transcripts. Without it Salmon drops them, they are missing from `quant.sf`, and their gene would get 0 silently.)
  - `pub fn salmon_quant(sample: &Sample, index_dir: &str, out_dir: &str, cfg: &RunConfig) -> Result<ToolCall, ConfigError>`. `outputs` list `<out_dir>/quant.sf` and `<out_dir>/aux_info/meta_info.json`.
  - `pub fn gene_counts_from_salmon(quant_sf: &str, tx2gene: &HashMap<String, String>) -> Result<BTreeMap<String, u64>, CountError>` (columns `Name`, `NumReads`).

- [ ] **Step 1: Add the tool to the environment.** Add `salmon = "==2.8.0"` to `[dependencies]` in `pixi.toml` (use `[workspace]`, never `[project]`). Run `pixi install` in the background with a log, as the hardening round did (`qRNAmap_data/pixi_install_harden.log`). Confirm the solve works on linux-64. If 2.8.0 does not solve, stop and report; do not pick another version. Then run `salmon index --help`, `salmon quant --help` and write the exact spelling of `-t -i -p`, `-l`, `-1 -2 -r`, `--fldMean`, `--fldSD`, `-o`, `--validateMappings` into the task report. Also test with the real binary and report: (a) whether a FASTQ path containing a comma is accepted (spec: Salmon rejects commas in paths; if it does, add `ConfigError::CommaInPath` to `salmon_quant` with a test, as HISAT2 does); (b) whether `-o` may already exist; (c) whether `--keepDuplicates` is accepted by `salmon index` in 2.8.0; (d) whether `.fastq.gz` input is read (Task 3 uses this). Run a tiny `salmon quant` on the `large.rs` fixture later (Task 3) to confirm the `meta_info.json` field names (`num_processed`, `num_mapped`, `percent_mapped`).
- [ ] **Step 2: Write failing tests in `crates/qrnamap-core/tests/toolcall.rs`:**
  - `salmon_quant_paired_unstranded`: args contain `-l IU`, then `-1 a.fq`, `-2 b.fq`, `-p 1`, `-o out`.
  - `salmon_quant_two_lanes_paired`: args hold `-1 l1_1.fq l2_1.fq` and `-2 l1_2.fq l2_2.fq` in that order (R1 list, then R2 list).
  - `salmon_quant_strandedness_paired`: Forward gives `ISF`, Reverse gives `ISR`.
  - `salmon_quant_single_end`: `-l U`, `-r a.fq`, `--fldMean 200 --fldSD 20` (f64 Display, as kallisto); Forward `SF`, Reverse `SR`.
  - `salmon_quant_single_end_needs_fragment_length`: `ConfigError::FragmentLengthRequired`. Also `InvalidFragmentLength` for NaN or 0.
  - `salmon_quant_rejects_bad_input`: `InvalidThreads`, `NoReads`, `MateCountMismatch`.
  - `salmon_index_args`: `index -t t.fa -i idx -p 4 --keepDuplicates`; threads 0 gives `InvalidThreads`.
  - `quantifier_name`: `Kallisto.name() == "kallisto"`, `Salmon.name() == "salmon"`.
- [ ] **Step 3: Write failing tests in `crates/qrnamap-core/tests/counts.rs`:**
  - `salmon_gene_sum_rounds_once`: two transcripts of one gene with `NumReads` 0.4 and 0.4 give 1 (sum 0.8); 0.25 + 0.25 gives 1 (0.5 rounds away from zero).
  - `salmon_unknown_transcript`, `salmon_missing_column` (`MissingColumn("NumReads")`), `salmon_bad_numreads` (negative, `nan`, `inf`, text give `BadRow`), `salmon_zero_gene_is_present_with_zero`.
  - `salmon_and_kallisto_agree_on_same_numbers`: same estimates in both table layouts give the same map.
- [ ] **Step 4: Run both test files. Expected: FAIL (items not defined).**
- [ ] **Step 5: Implement.** Export `salmon_index`, `salmon_quant` and `gene_counts_from_salmon` from `crates/qrnamap-core/src/lib.rs` (the `pub use` list at about lines 16-23). Add the enum variant and `name()` in `config.rs` (fix the doc comment). In `toolcall.rs` add the two `Tool` variants, `program()` arms, the two builders; reuse `check_threads`, `check_reads` and the fragment-length validation (factor the validation into a private helper shared with `kallisto_quant` without changing kallisto output). In `counts.rs` factor the existing sum-and-round body of `gene_counts_from_abundance` into a private helper taking the id and count column names, and call it from both public functions. Kallisto tests must stay green and unchanged.
- [ ] **Step 6: Run `cargo test -p qrnamap-core`. Expected: all PASS. Then run the wasm32 build check for `qrnamap-core` the way the ledger records it (hardening ledger, wasm check). Expected: builds.**
- [ ] **Step 7: Commit** (`feat(core): Salmon builders, library types and gene sum`).

---

### Task 2: Summary format gets a `quantifier` column

**Files:**
- Modify: `crates/qrnamap-cli/src/pipeline.rs` (`SUMMARY_HEADER`, `write_sample_summary`, `read_sample_summary`, `build_mapping_summary`, strand-check header at the `strand_check_suggestion` line), callers in `quant_sample` and `align_sample`.
- Test: `crates/qrnamap-cli/tests/{mapping_summary,cli,run}.rs` and every test that asserts the 7-column header or calls `write_sample_summary`. Find them with `grep -rn "write_sample_summary\|qrnamap_version\|f\[[0-9]\]" crates/qrnamap-cli/tests`. Known sites: `tests/matrix.rs:54`, `tests/run.rs:932`, `tests/mapping_summary.rs` (cell indexes such as `f[4]` shift by one).

**Interfaces:**
- Produces: `SUMMARY_HEADER = "sample\tmode\tquantifier\treads_processed\tfraction_mapped\tfraction_assigned\tstrandedness\tqrnamap_version"` (8 columns). `write_sample_summary(out_dir, sample, mode, quantifier: Option<Quantifier>, summary, strandedness)`. The cell is `kallisto`, `salmon`, or `none` (align).

- [ ] **Step 1: Update the existing tests first** to the 8-column layout (column 3 `kallisto` for quant, `none` for align). Add `old_seven_column_summary_is_rejected`: a hand-written old file makes `build_mapping_summary` return `BadTable` whose detail names the header. Add `mixed_quantifier_rows_keep_their_tool`: one kallisto row and one salmon row in a built `mapping_summary.tsv.zst` (write the per-sample files directly with `write_sample_summary`).
- [ ] **Step 2: Run them. Expected: FAIL.**
- [ ] **Step 3: Implement.** Change the header, writer, reader (`row.len() != 8`), every place that indexes summary cells (`build_mapping_summary` and the strand-check suggestion column; find them with `grep -n "row\[" crates/qrnamap-cli/src/pipeline.rs`), and both sample functions (kallisto passes `Some(Quantifier::Kallisto)` for now; align passes `None`).
- [ ] **Step 4: Run `cargo test -p qrnamap-cli`. Expected: PASS, with the tool-gated tests run under `QRNAMAP_REQUIRE_TOOLS=1` and the pixi PATH.**
- [ ] **Step 5: Commit** (`feat(cli): quantifier column in the sample summary`).

---

### Task 3: Salmon in `prepare_quant` and `quant_sample`

**Files:**
- Modify: `crates/qrnamap-cli/src/pipeline.rs`
- Test: `crates/qrnamap-cli/tests/pipeline_quant.rs`

**Interfaces:**
- Consumes: Task 1 builders and `gene_counts_from_salmon`; Task 2 summary writer.
- Produces:
  - `QuantContext` has no new field: use `ctx.config.quantifier`. `index` is a file (kallisto, `transcripts.kidx`) or a directory (Salmon, `<out>/salmon_index`).
  - `prepare_quant` builds only the index of `config.quantifier`.
  - `PipelineError::PairCountMismatch { sample, input, processed, tool: &'static str }` (rename the field, name the tool in the message; update the match that uses `kallisto: 4` at about pipeline.rs:1442). `NoReadsPseudoaligned` gets a `tool` field and names it. `check_processed(sample, input, processed, tool)`.
  - `fn parse_salmon_meta(text: &str, path: &Path) -> Result<(u64, u64, f64), PipelineError>`: `(num_processed, num_mapped, percent_mapped)`. Private is fine. It rejects a key that occurs twice and a missing key (`MissingField`).
  - `remove_old_outputs` also removes `quant.sf`, `quant.sf.partial` if any, `aux_info/meta_info.json`, `cmd_info.json`, and the Salmon per-sample subdirectories it created (`aux_info`, `logs`, `libParams`), `lib_format_counts.json` and `libParams/flenDist.txt`; kallisto files stay in the list.

- [ ] **Step 1: Write failing tests** (tool-gated ones use the existing skip helper in `pipeline_quant.rs`):
  - `salmon_quant_sample_matches_truth_on_fixture` (gated): prepare with Salmon, quantify the small fixture, assert gene counts within the tolerance the kallisto test uses, summary row says `salmon`, `reads_processed` equals the input pairs.
  - `salmon_zero_mapped_reads_is_named_error` (gated): reads that match nothing give `NoReadsPseudoaligned` with the tool `salmon`. Salmon 2.x writes its output and then exits non-zero when no fragment maps (salmon MIGRATION.md), so the error comes from the `ExecError::Failed` path, not from the parser.
  - `salmon_pair_count_mismatch_names_tool` (unit, in a `#[cfg(test)] mod tests` inside `pipeline.rs`, because the function is private): `check_processed("s", 10, 9, "salmon")` message contains `salmon`.
  - `parse_salmon_meta_ok_and_errors` (unit, same `mod tests`): good JSON; missing key; duplicated key (a nested object repeating `num_processed`); non-numeric value.
  - `switching_tool_removes_old_outputs` (unit, no tools): create fake `abundance.tsv`, `run_info.json`, `quant.sf`, `aux_info/meta_info.json`, summary and genecounts in a sample dir; `remove_old_outputs` removes all of them.
  - `kallisto_gene_counts_unchanged` (gated, regression): quantify the large fixture with kallisto and compare with the counts stored in a new golden file `crates/qrnamap-cli/tests/golden/large_kallisto_counts.tsv`. Create the golden file only with a binary built from commit `1300552` (the plan commit, before Task 1) in a separate checkout (`git worktree add <scratch dir> 1300552`). Do not build it on this branch: Task 1 already refactored the shared kallisto code. Record the commit, the kallisto build string (`kallisto version` and the conda build from `pixi.lock`) and the command in the report. If the kallisto build changed with the Task 1 `pixi.lock` update, say so.
  - `check_samples_dry_run_builds_salmon_call` (unit): with Salmon config and a sample whose path contains a comma (only if Task 1 found that Salmon rejects commas), `check_samples` returns the error from the Salmon builder, while the same sample with kallisto config passes. If Task 1 found no Salmon-only check, delete this test and keep the branch covered by the gated tests.
- [ ] **Step 2: Run them. Expected: FAIL.**
- [ ] **Step 3: Implement.** `check_samples`: branch on `cfg.quantifier` (dry call to `kallisto_quant` or `salmon_quant` with placeholders). `prepare_quant`: for Salmon, `plain_copy` the transcripts as today, call `salmon_index(.., &salmon_index_dir, &config)`; run it with `run_call`; remove a stale `salmon_index` directory first (R24 order: check inputs, remove, build). `quant_sample`: same staging and validation as now (pass `gzip_ok` per the Step 1 finding: if Salmon 2.8.0 reads `.gz`, keep `true`; otherwise `false` so the layer decompresses). Run `salmon_quant`; on `ExecError::Failed`, if `<sample>/aux_info/meta_info.json` exists and parses with `num_mapped == 0`, return `NoReadsPseudoaligned { tool: "salmon" }`, otherwise `Exec` (R24 removed any old `meta_info.json` first, so a stale file cannot match); if the run succeeds with `num_mapped == 0`, return the same error; parse `quant.sf` and `aux_info/meta_info.json`; `check_processed` against the input pairs; `fraction_mapped = fraction_assigned = percent_mapped / 100`; write gene counts and the summary with `Some(Quantifier::Salmon)`. Kallisto branch keeps its code, only passing `tool: "kallisto"` to the renamed errors.
- [ ] **Step 4: Run `cargo test -p qrnamap-cli` with `QRNAMAP_REQUIRE_TOOLS=1`. Expected: PASS, and the regression test shows identical kallisto counts.**
- [ ] **Step 5: Commit** (`feat(cli): salmon prepare and quant`).

---

### Task 4: Flag, marker, `--task-id` and SLURM

**Files:**
- Modify: `crates/qrnamap-cli/src/main.rs` (flag on the `Quant` and `Run` variants only, not on `RunFlags`, because `RunFlags` is also flattened into `Align` at about line 105; set `cfg.quantifier` after `run.config()`; help text lines 15, 67, 303), `crates/qrnamap-cli/src/run.rs` (`MARKER_HEADER`, `write_prepared_marker`, `parse_marker`, `check_prepared`, `prepared_quant`, `require_file`, `required_tools`, `slurm_plan`, `remove_old_run_outputs`), `crates/qrnamap-cli/src/slurm.rs` if it names kallisto.
- Test: `crates/qrnamap-cli/tests/{run,slurm,cli}.rs`

**Interfaces:**
- Produces: `--quantifier <kallisto|salmon>` (clap `value_enum`, default `kallisto`). `required_tools(mode: Mode, quantifier: Quantifier, with_de: bool)`. Marker file: header `qrnamap-prepared\t2`, new line `quantifier\t<kallisto|salmon|none>` after `mode`. `write_prepared_marker(out, mode, quantifier: Option<Quantifier>, inputs)` and `check_prepared(out, mode, quantifier: Option<Quantifier>, inputs)`; `None` is for align and writes `none`. Update the callers (run.rs about lines 431, 448, 508, 553) and `tests/run.rs:595`, which asserts `qrnamap-prepared\t1`.

- [ ] **Step 1: Write failing tests:**
  - `marker_records_quantifier_and_refuses_other` (tests/run.rs): prepare with kallisto, then `check_prepared` with Salmon returns `NotPrepared` whose reason contains both `kallisto` and `salmon`.
  - `marker_version_1_is_refused_with_advice`: a hand-written version-1 marker gives a reason that says to prepare again.
  - `marker_for_align_writes_none`.
  - `task_id_with_other_quantifier_refuses` (CLI-level, may use the existing fake/dry harness in `tests/run.rs`): `--step final` and `--task-id` on a kallisto-prepared directory with `--quantifier salmon` exit non-zero.
  - `quantifier_with_align_mode_is_error` (tests/cli.rs): `qrnamap run --mode align --quantifier salmon` exits with a message that names both options. The default value must not trigger this: make the flag `Option<QuantArg>` in clap and resolve the default after the mode check.
  - `slurm_forwards_quantifier` (tests/slurm.rs): the generated prepare, samples and final scripts contain `--quantifier salmon`; `required_tools` with Salmon lists `salmon` and not `kallisto`; with kallisto the reverse. Update existing `required_tools` callers.
  - `slurm_rate_note_names_kallisto_only`: the `slurm.rs` text that prints per-million-pair rates (about lines 29-30 and 272) says the rates were measured with kallisto; for Salmon it says they are not measured. `prepared_quant_accepts_a_directory_index`: Salmon prepared directory with `salmon_index/` passes `require_file`-style check; a missing directory gives `NotPrepared`.
- [ ] **Step 2: Run them. Expected: FAIL.**
- [ ] **Step 3: Implement** as in the Interfaces block. `prepared_quant` takes the quantifier from `config.quantifier`, checks the marker with it, and checks `transcripts.kidx` (file) or `salmon_index` (directory; use a `require_dir` twin of `require_file`). `RunFlags::config` sets `quantifier`. `run` errors when `--quantifier` is given and `--mode align`; `qrnamap align` has no such flag. In `slurm_plan`, forward `--quantifier` in `run_args` and `task_args` only when the mode is quant. `remove_old_run_outputs` needs no change beyond Task 3's removal list.
- [ ] **Step 4: Run `cargo test --workspace` with `QRNAMAP_REQUIRE_TOOLS=1` and the pixi PATH. Expected: PASS. Run the wasm build command from `pixi.toml` directly (R11; CI runs only `pixi run test` and `pixi run wasm-check`). Expected: builds. `cargo clippy` is optional and not a CI gate.**
- [ ] **Step 5: Commit** (`feat(cli): --quantifier flag, marker and SLURM support`).

---

### Task 5: Tool-gated end-to-end tests on the large fixture

**Files:**
- Modify: `crates/qrnamap-cli/tests/large.rs` (`QUANT_TOOLS` is a kallisto-only const at about line 110), `crates/qrnamap-cli/tests/run.rs`.

- [ ] **Step 1: Write failing gated tests:**
  - `salmon_recovers_counts_on_large_fixture`: counts within 10 percent of truth for genes above a stated count floor (use the floor the kallisto test uses).
  - `salmon_two_lanes_equal_concatenated`: counts from two lanes equal counts from the concatenated file (exact equality; if Salmon is not deterministic at the single-thread setting, report the measured difference and assert within 1 count; do not loosen without reporting).
  - `salmon_forward_and_reverse_stranded`: a stranded fixture run gives the expected sense counts with the right `-l`, and a wrong strandedness gives either clearly lower assigned counts or the `NoReadsPseudoaligned` error, as the kallisto test at `large.rs:237-261` accepts (report the numbers).
  - `salmon_single_end_run`: single-end sample with `--fragment-length 200,20` finishes and recovers counts.
  - `run_end_to_end_with_salmon_quantifier`: `qrnamap run --quantifier salmon` produces `counts.tsv.zst`, `genes.tsv.zst`, `mapping_summary.tsv.zst` with `salmon` in the quantifier column.
- [ ] **Step 2: Run them under `QRNAMAP_REQUIRE_TOOLS=1`. Expected: FAIL first (if a test cannot fail first because the code exists, state that in the report and show it passing).**
- [ ] **Step 3: Fix any real defect they expose in the code (not in the tests).** Add a separate `SALMON_TOOLS` list for the Salmon tests. Do not add `salmon` to the shared `QUANT_TOOLS` (large.rs:110), or kallisto tests would skip, or fail under `QRNAMAP_REQUIRE_TOOLS=1`, when only salmon is missing.
- [ ] **Step 4: Run the full workspace tests. Expected: PASS.**
- [ ] **Step 5: Commit** (`test: Salmon end-to-end tests on the large fixture`).

---

### Task 6: Documentation and AMB validation

**Files:**
- Modify: `README.md` (kallisto mentions near lines 4, 61, 155, 168, 172, 191-201, 233 ("No Salmon"), 275-290; document `--quantifier`, the index directory, the 8-column summary, the marker version, Salmon version 2.8.0), `validation/README.md`, new `validation/amb/salmon.md`, `docs/superpowers/specs/2026-11-19-qrnamap-design.md` (note in the Salmon lines that the addendum supersedes them), and the main `--help` text.

- [ ] **Step 1: Update the docs** as above. State limits plainly: transcripts-only index, no decoys, no bootstrap, no align-mode Salmon, summary format changed.
- [ ] **Step 2: AMB run.** Using the pixi PATH and the AMB data under `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb`, run `qrnamap run --quantifier salmon` as a SLURM job set (`run --slurm`, default `--mem-gb 8`, no `--tmpdir` unless a measured reason exists, so jobs use node-local `$SCRATCH`). Compare Salmon gene counts with the kallisto counts from the earlier AMB run and with the authors' GSE80357 tables, with the same script and metrics as the GEO comparison (`validation/geo/`). Write measured numbers only (correlation, number of genes, DE agreement if the GEO script reports it) in `validation/amb/salmon.md`, with the job IDs, wall time and MaxRSS from `sacct`. If a number is not measured, write that it is not measured.
- [ ] **Step 3: Run the full workspace tests, clippy, fmt and the wasm check once more. Expected: all pass.**
- [ ] **Step 4: Commit** (`docs: Salmon support and AMB comparison`).
