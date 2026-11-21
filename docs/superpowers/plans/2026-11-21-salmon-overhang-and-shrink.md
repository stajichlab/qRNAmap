# Salmon overhang options, shrinkage packages and help text: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add opt-in Salmon overhang options, install `apeglm` and `ashr`, fix the `qrnamap de --help` text, and measure the Salmon options on the six AMB samples.

**Architecture:** Two typed fields in `RunConfig`, one builder change in `salmon_quant`, two flags on `quant` and `run`, forwarded by `run --slurm`. Pixi gets two R packages. A measurement task reruns AMB.

**Tech Stack:** Rust workspace, Salmon 2.8.0, R 4.5 + DESeq2 1.50.2, pixi.

**Spec:** `docs/superpowers/specs/2026-11-21-salmon-overhang-and-shrink-design.md` (read it first). Related: `docs/superpowers/specs/2026-11-20-salmon-support-design.md`, `validation/amb/salmon.md`.

## Global Constraints

- Defaults unchanged: with neither new flag, the Salmon command line and all outputs are identical to `main` (fe1ce5e). Kallisto is untouched.
- The DESeq2 pin `bioconductor-deseq2 ==1.50.2` and `r-base` do not change.
- Run tools with `PATH=<worktree>/.pixi/envs/default/bin:$PATH`. Never `pixi run` (R11). `CARGO_TARGET_DIR=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-overhang`.
- The new worktree has no pixi env yet: Task 1 creates it with `pixi install` in the background (log under `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/`), as earlier rounds did.
- Pure crates build for wasm32. No JSON dependency in core. R5, R18, R24, R27 as before.
- Flag names: `--salmon-softclip-overhangs` (switch), `--salmon-min-score-fraction F` with `F` finite and in (0, 1]. Salmon spellings: `--softclipOverhangs`, `--minScoreFraction`.
- SLURM jobs use default `--mem-gb 8` and node-local `$SCRATCH`; no `--tmpdir` unless a measured reason exists.
- Commit trailers on every commit: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01CoSBcPmuXVDpBrrPNsnYtc`.
- Docs style: Simplified Technical English. Measured numbers only. No `cargo fmt` surprises: the repo is rustfmt-clean now, so run `cargo fmt --all -- --check` before each commit and keep it clean.

## Review Focus

- A flag combination that is nonsense: `--salmon-min-score-fraction 0`, `1.5`, `NaN`, `-0.1` (must be an error naming the option, before any tool runs).
- The flags with `--quantifier kallisto` or the default, and with `--mode align` (error naming both options; the default must not trigger it).
- `run --slurm` with the flags: the generated samples script carries them; without the flags the scripts are byte-identical to `main`.
- Salmon default command line unchanged when neither flag is set (test the exact argv).
- Adding the R packages must not change the DESeq2 version printed by `qrnamap de`.

## File Structure

- Modify `pixi.toml`, `pixi.lock`.
- Modify `crates/qrnamap-core/src/config.rs`, `toolcall.rs`, tests `crates/qrnamap-core/tests/toolcall.rs`.
- Modify `crates/qrnamap-cli/src/main.rs`, `run.rs` (slurm plan), tests `tests/cli.rs`, `tests/run.rs`, `tests/slurm.rs`, `tests/large.rs`, `tests/de_run.rs`.
- Modify `README.md`, `validation/amb/salmon.md`; add `validation/amb/run_amb_salmon_overhang.sh` and results under `validation/amb/results_salmon_overhang/`.

---

### Task 1: Shrinkage packages and help text (item 4)

**Files:** `pixi.toml`, `pixi.lock`, `crates/qrnamap-cli/src/main.rs` (the `De` variant doc comments, `--samples`, `--shrink`), `README.md`, `crates/qrnamap-cli/tests/de_run.rs`.

**Produces:** a working pixi env in this worktree with `apeglm` and `ashr`; help text; gated tests.

- [ ] **Step 1: Add `bioconductor-apeglm` and `r-ashr` to `[dependencies]`.** Run `pixi install` in the background with a log. Confirm the solve on linux-64. Confirm `bioconductor-deseq2` stays at `1.50.2` and `r-base` stays on 4.5.x (compare `pixi.lock` before and after for those packages; report any other package that moves). If a package does not solve for R 4.5, stop and report BLOCKED with the solver output; do not change R or DESeq2.
- [ ] **Step 2: Write failing tests in `tests/de_run.rs`:** `shrink_apeglm_writes_shrunk_results` and `shrink_ashr_writes_shrunk_results` (gated on the R package being installed: skip without it, fail under `QRNAMAP_REQUIRE_TOOLS=1`; follow how the existing `--shrink normal` test gates and checks `PREFIX.results.shrunk.csv`). Assert the file exists, has the expected columns for each method and that `log2FoldChange` values are finite for most rows. Also `de_version_line_unchanged_by_shrink_packages`: `qrnamap de` still prints `DESeq2 1.50.2`.
- [ ] **Step 3: Run them. Expected: FAIL or SKIP before the install; PASS after Step 1 if the code already supports these methods (state in the report which).** If the generated R script needs a change for apeglm or ashr to run (for example a coefficient name), fix it in `crates/qrnamap-core/src/de.rs` with a test in `crates/qrnamap-core/tests/de_script.rs`.
- [ ] **Step 4: Run `--shrink apeglm` and `--shrink ashr` on the AMB counts** (counts and sample sheet under `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb`, use the earlier AMB DE contrast from `validation/amb/`). Report the commands, the number of rows, and the correlation between normal, apeglm and ashr shrunk log fold changes (measured). Output stays out of git except a small summary table under `validation/amb/results_shrink/`.
- [ ] **Step 5: Help text.** `--samples`: "Sample sheet (TSV). Its columns hold the design variables." (improve using the sheet rules in README). `--shrink`: state that the pixi environment has `apeglm` and `ashr`, that `normal` needs no extra package, and the output file name. Update the README `de` section the same way. Remove "which pixi.toml does not install; only normal was run".
- [ ] **Step 6: Run `cargo test -p qrnamap-cli --test de_run` and `cargo test -p qrnamap-core` with `QRNAMAP_REQUIRE_TOOLS=1`. Expected: PASS.** Commit (`feat: apeglm and ashr in the pixi environment, de help text`).

---

### Task 2: Core: config fields and Salmon builder options

**Files:** `crates/qrnamap-core/src/config.rs`, `toolcall.rs`, `lib.rs` if needed; test `crates/qrnamap-core/tests/toolcall.rs`.

**Interfaces — Produces:** `RunConfig.salmon_softclip_overhangs: bool` (default `false`), `RunConfig.salmon_min_score_fraction: Option<f64>` (default `None`); `ConfigError::InvalidMinScoreFraction { value: f64 }` with a Display text that names `--salmon-min-score-fraction`; `salmon_quant` appends `--softclipOverhangs` when the bool is set and `--minScoreFraction <F>` when the option is set, before `-p`.

- [ ] **Step 1: Write failing tests:** `salmon_quant_default_argv_unchanged` (exact full argv for a paired unstranded sample equals the argv `main` produces; copy it from the existing test expectations), `salmon_quant_softclip_overhangs_flag`, `salmon_quant_min_score_fraction_flag` (`0.5` prints `0.5`; `1` prints `1`), `salmon_quant_both_flags_order` (both flags appear before `-p`), `salmon_quant_min_score_fraction_rejects_bad_values` (`0.0`, `-0.1`, `1.5`, `NaN`, `inf` give `InvalidMinScoreFraction`), `kallisto_quant_ignores_salmon_options` (the kallisto argv is identical with and without the fields set).
- [ ] **Step 2: Run `cargo test -p qrnamap-core --test toolcall`. Expected: FAIL (fields undefined).**
- [ ] **Step 3: Implement.** Update every `RunConfig { .. }` literal in the workspace that lists all fields (grep). Keep lines under 100 chars.
- [ ] **Step 4: Run `cargo test -p qrnamap-core`, `cargo fmt --all -- --check`, and the wasm32 build of `qseq-io`, `genome-catalog`, `qrnamap-core`. Expected: PASS and builds.** Commit (`feat(core): Salmon overhang and min-score-fraction options`).

---

### Task 3: CLI flags, SLURM forwarding and end-to-end tests

**Files:** `crates/qrnamap-cli/src/main.rs` (flags on `Quant` and `Run`, resolution, errors), `run.rs` (`slurm_plan`: `run_args` and `task_args`), tests `tests/cli.rs`, `tests/run.rs`, `tests/slurm.rs`, `tests/large.rs`.

**Interfaces — Consumes:** Task 2 config fields. **Produces:** flags `--salmon-softclip-overhangs` and `--salmon-min-score-fraction <F>` on `quant` and `run`.

- [ ] **Step 1: Write failing tests:**
  - `salmon_flags_rejected_with_kallisto` and `..._with_default_quantifier`: error names the flag and `--quantifier salmon`.
  - `salmon_flags_rejected_with_align_mode` (`run --mode align`): error names both options; `run --mode align` without the flags still passes the check.
  - `salmon_min_score_fraction_invalid_value_is_error`: `0`, `1.5`, `abc` fail with a message naming the option, before any file is written (R20).
  - `slurm_forwards_salmon_flags`: the samples script and the run args contain `--salmon-softclip-overhangs` and `--salmon-min-score-fraction 0.5`; `slurm_without_flags_is_unchanged`: scripts for a Salmon plan without the flags equal the scripts `main` generates (compare against a stored expected string taken from `main`; document how).
  - Gated end-to-end on the large fixture (`tests/large.rs`, use `SALMON_TOOLS`): `salmon_overhang_flags_run_and_record_command`: quant with both flags finishes, recovers counts within the same tolerance as the default Salmon test, and `<out>/<sample>/cmd_info.json` contains `softclipOverhangs` and `minScoreFraction`; `salmon_defaults_have_no_overhang_options_in_cmd_info`: without flags `cmd_info.json` contains neither (check the file's actual field names with a real run and report them).
  - Measured: Salmon 2.8.0 records neither flag in `cmd_info.json`. The assertion was replaced by `salmon_overhang_flags_accepted_counts_within_tol` (ledger ruling, commit 239ae3d).
- [ ] **Step 2: Run them. Expected: FAIL.**
- [ ] **Step 3: Implement** in the same pattern as `--quantifier` (flags on `Quant` and `Run` only, resolved after the mode check). Update help text and README flag list. Forward in `slurm_plan` only when set.
- [ ] **Step 4: Run `cargo test --workspace` with `QRNAMAP_REQUIRE_TOOLS=1`, `cargo fmt --all -- --check`. Expected: PASS.** Commit (`feat(cli): --salmon-softclip-overhangs and --salmon-min-score-fraction`).

---

### Task 4: AMB measurement and docs

**Files:** `validation/amb/run_amb_salmon_overhang.sh`, `validation/amb/results_salmon_overhang/`, `validation/amb/salmon.md`, `README.md`.

- [ ] **Step 1: Run the six AMB samples with Salmon** for two new settings: (A) `--salmon-softclip-overhangs`, (B) `--salmon-softclip-overhangs --salmon-min-score-fraction 0.5`. Reuse the infrastructure of `validation/amb/run_amb_salmon.sh` and `collect_salmon.py` (read them first). Build the binary from this branch. Use `qrnamap run --slurm` (default `--mem-gb 8`, no `--tmpdir`, `REPO` passed explicitly, never `BASH_SOURCE`). Outputs go under `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb/`. The existing default-setting Salmon run and the kallisto run are the baselines.
- [ ] **Step 2: Compare** with the existing scripts (`zero_genes_by_length.py`, the compare scripts used for `results_salmon`): for each setting, per sample and pooled: mapped fraction; number of zero-count genes; genes zero in Salmon and non-zero in kallisto, and how many of those have all transcripts shorter than the reads; RPL41A/B and YOR302W counts; gene-count Spearman with kallisto and with the authors' GSE80357 table; DE agreement with the authors' DE (Jaccard) as the earlier comparison reports. Commit small TSVs.
- [ ] **Step 3: Write `validation/amb/salmon.md` section "Overhang options"** with only measured numbers, job IDs, `sacct` wall time and MaxRSS, commands. State plainly what improved and what did not. If a setting recovers some genes but changes agreement with the authors' tables, give both numbers. No recommendation to change the default unless the numbers show it; list it as the user's decision. Update the README Salmon limit text and the "OPEN USER DECISION" note to point to the new numbers.
- [ ] **Step 4: Run the full workspace tests, `cargo fmt --all -- --check`, and the wasm check. Expected: all pass.** Commit (`docs: Salmon overhang options measured on AMB`).
