# Salmon overhang options, shrinkage packages and help text: design

Date: 2026-11-21. Status: approved by the user on 2026-11-21 (design given in conversation).
Parent: `2026-11-20-salmon-support-design.md` and `validation/amb/salmon.md`.

## Intent

Salmon with default options gives 0 reads to transcripts shorter than the read length. On the AMB data (101 nt reads), 288 to 306 genes per sample are 0 in Salmon and non-zero in kallisto. RPL41A and RPL41B are among them. A 500,000-pair subset of one sample showed that `--softclipOverhangs` recovers part of these reads. The user wants these options available and wants the effect measured on the full data.
The user also deferred two small items: the `--samples` text in `qrnamap de --help`, and the packages `apeglm` and `ashr` so that `--shrink apeglm|ashr` can run.

## Goals

1. `--salmon-softclip-overhangs` (switch) and `--salmon-min-score-fraction F` on `qrnamap quant` and `qrnamap run`. Both are off by default. With neither flag the Salmon command and the output stay as they are now.
2. `run --slurm` forwards the flags to the generated scripts.
3. A full re-run of the six AMB samples with: defaults (existing result), `--salmon-softclip-overhangs`, and `--salmon-softclip-overhangs --salmon-min-score-fraction 0.5`. Report zero-count genes, agreement with kallisto, agreement with the authors' GSE80357 tables, and the RPL41A/B counts. Measured numbers only.
4. `qrnamap de --help` describes `--samples` and no longer says that `apeglm` and `ashr` are not installed, if the packages work.
5. `bioconductor-apeglm` and `r-ashr` in `pixi.toml`. `--shrink apeglm` and `--shrink ashr` run on the AMB counts and give output.

## Non-goals

- Changing a default. The default changes only if the user decides so after reading the measured numbers.
- A generic `--salmon-arg` passthrough. The explicit flags stay validated and argv-only (R18).
- Genome decoys, bootstrap, align-mode Salmon, the browser build, `apeglm`/`ashr` in the browser.
- Other Salmon tuning options (`--softclip`, `--validateMappings`).

## Anti-goals

- Any change to kallisto behavior or to Salmon output when neither flag is set.
- A change of the DESeq2 version when the shrinkage packages are added.

## Approaches

1. **Explicit flags (chosen).** Two typed fields in `RunConfig`, validated, forwarded by SLURM. Cost: one more field per option if more options come.
2. Generic passthrough `--salmon-arg X`. Less code, no validation, settings invisible in help. Rejected.
3. Config file for tool options. More machinery than two options need. Rejected.

## Design

- `qrnamap-core::config::RunConfig`: `salmon_softclip_overhangs: bool` (default false), `salmon_min_score_fraction: Option<f64>` (default None).
- `qrnamap-core::toolcall::salmon_quant`: append `--softclipOverhangs` when set, and `--minScoreFraction <F>` when set, before `-p`. `F` must be finite and in (0, 1]. Otherwise `ConfigError::InvalidMinScoreFraction`. f64 formatting as the fragment-length values (`0.5` prints `0.5`). The flags apply to quant only. The index command is unchanged.
- CLI: the flags are on the `Quant` and `Run` variants only. A flag with `--quantifier kallisto` (or the default) or with align mode is an error that names both options. The default resolution follows the `--quantifier` pattern (`Option`/`bool` resolved after the mode check).
- SLURM: `slurm_plan` puts the flags in `run_args` and `task_args` when set and the quantifier is Salmon. Prepare does not need them but receives them with `run_args` (as `--quantifier` does).
- The `.prepared` marker does not record the flags: they do not change the index. Salmon 2.8.0 does not write either flag to any output file (measured: `cmd_info.json` has only salmon_version, index, libType, output, mates1, mates2, unmatedReads, threads, sketch). The generated SLURM scripts and the qrnamap command line are the record of the options.
- Pixi: add `bioconductor-apeglm` and `r-ashr` (conda-forge or bioconda, R 4.5). Check the solve. The DESeq2 pin `==1.50.2` must not change. Record the versions in `qrnamap de` version output if it lists packages (builder checks).
- Help text: `--samples` gets "Sample sheet (TSV); its columns hold the design variables." (builder improves the wording from the sheet rules). `--shrink` text states which packages the pixi environment has and that `normal` needs none.

## Constraints

- Pure crates build for wasm32. No JSON dependency in core.
- Never `pixi run`. Use the env bin on PATH.
- R5 (tool tests skip or fail under `QRNAMAP_REQUIRE_TOOLS=1`), R18, R24, R27 as before.
- Jobs use the default `--mem-gb 8` and node-local `$SCRATCH` (no `--tmpdir` unless a measured reason exists).

## Left to the builder

- Exact solver choices for `apeglm` and `ashr` builds, and whether R 4.5 packages exist for linux-64.
- The test design for `--shrink apeglm|ashr` (gated on the R package being present).
- The exact wording of the help text.
