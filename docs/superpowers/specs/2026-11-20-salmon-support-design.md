# Salmon quantifier support: design addendum

Date: 2026-11-20. Status: draft for approval.
Parent design: `2026-11-19-qrnamap-design.md` (sections 157, 158, 204, 207 mention Salmon).

## Intent

The user asked for Salmon as a second quantifier next to kallisto. Kallisto stays the default. The user wants to compare the two tools on the same data and to use Salmon where it fits better.

## Goals

1. `--quantifier kallisto|salmon` on `qrnamap quant` and `qrnamap run`. The default is `kallisto`.
2. Salmon output enters the same downstream files as kallisto output: `counts.tsv.zst`, `genes.tsv.zst`, `<sample>.summary.tsv.zst`, `mapping_summary.tsv.zst`. The DE step needs no change.
3. Salmon works in the SLURM scripts (`run --slurm`) with the same safety rules as kallisto (R24, R27, R28).
4. Paired-end, single-end, multi-lane samples, and all three strandedness settings work.

## Non-goals

- Genome decoys in the Salmon index. The user chose option A: transcripts only. Decoys are a later step.
- Salmon in the browser (WASM) build.
- Bootstrap or Gibbs sampling, `--validateMappings` tuning, GC or sequence bias correction. Use Salmon defaults and note the tool version.
- Salmon in align mode. Align mode keeps HISAT2 and featureCounts.

## Anti-goals

- Any change to kallisto behavior or output. A regression test must show identical kallisto gene counts before and after.
- A silent default change. Existing commands and scripts must give the same result as today.

## Approaches

1. **Branch on `Quantifier` inside the existing code (chosen).** Add Salmon argv builders, a `quant.sf` gene-sum parser and a summary parser as small pure functions in `qrnamap-core`. The CLI picks between the two tools in `prepare_quant`, `quant_sample` and the SLURM tool list. Kallisto code stays as it is. Cost: `pipeline.rs` grows and gets branches.
2. Refactor both tools behind a backend trait first. Cleaner for a third tool and for the browser. It touches validated kallisto code. Deferred until the browser plan or a third tool needs it.
3. Copy the kallisto path into a separate Salmon pipeline. About 300 duplicated lines. Every later fix is needed twice. Rejected.

The user chose option 1 by default ("I will go with my pick unless you object").

## Design

### Core crate (`qrnamap-core`, pure, must build for wasm32)

- `config.rs`: add `Quantifier::Salmon`. Fix the doc comment.
- `toolcall.rs`:
  - `Tool::SalmonIndex`, `Tool::SalmonQuant`.
  - `salmon_index(transcripts, index_dir, cfg) -> ToolCall`: `salmon index -t <transcripts.fa> -i <index_dir> -p <threads>`.
  - `salmon_quant(sample, index_dir, out_dir, cfg) -> Result<ToolCall, ConfigError>`.
    - Paired: `-l <lib> -1 <all R1 files> -2 <all R2 files>`. Salmon takes lanes as a list after each flag. This differs from kallisto, which takes adjacent pairs.
    - Single-end: `-l <lib> -r <files>` plus `--fldMean <mean> --fldSD <sd>`. The existing fragment length option supplies the values. Same validation as kallisto (finite, > 0, required).
    - Library type is explicit, never `-l A`. Paired: `IU` (unstranded), `ISF` (forward), `ISR` (reverse). Single-end: `U`, `SF`, `SR`. The `Strandedness` enum maps to these.
    - `-p <threads>`. Output directory `-o <out_dir>`.
  - Same `ConfigError` checks as kallisto: threads >= 1, mate counts equal, no commas in paths.
- `counts.rs`:
  - `gene_counts_from_salmon(quant_sf_text, tx2gene) -> ...`: sum `NumReads` per gene through tx2gene, round once, half away from zero. Same output type and rounding as `gene_counts_from_abundance` (kallisto). Share the summing and rounding code where possible.
  - `parse_salmon_meta(meta_info_json) -> SalmonMeta { num_processed, num_mapped, percent_mapped }`. The reader must fail with a clear error if a field is missing. Builder: confirm the field names against real Salmon 2.8.0 output.
  - Correction after builder check: `parse_salmon_meta` lives in the CLI crate beside the existing `json_number` hand parser (pipeline.rs), not in core. Salmon in the browser is a non-goal, and core has no JSON dependency. The parser must check that each key occurs once, because `meta_info.json` has nested objects.

### CLI crate (`qrnamap-cli`)

- `--quantifier` on `quant` and `run`. Parsed to `Quantifier`. Default kallisto.
- Index directory: `<out>/salmon_index/` (kallisto keeps `<out>/transcripts.kidx`). `prepare_quant` builds the one that matches the quantifier. It does not build both.
- `quant_sample` for Salmon:
  - Validate the input like the kallisto path does. Read every lane to the end, check mate counts, count input pairs.
  - Run `salmon quant`. Read `<out>/<sample>/quant.sf` and `aux_info/meta_info.json`.
  - Compare `num_processed` with the input pair count. A mismatch fails the sample (same rule as the kallisto C1 check).
  - Write per-sample gene counts in the same form as kallisto so `assemble` is unchanged.
- Gzip and zstd input: the builder checks whether Salmon reads `.gz` directly. Zstd is staged by our own layer exactly as kallisto input is today. If Salmon cannot read gz in 2.8.0, stage it like kallisto.
- Stale output (R24): remove old Salmon output for the sample first, in the same order as kallisto (check inputs, remove, build).
- `.prepared` marker (R27): add the quantifier to the stamped inputs. `--task-id` and `--step final` refuse a marker that names a different quantifier.
- Summary files: `SUMMARY_HEADER` gets a `quantifier` column. This changes the file format for kallisto too. The builder updates the readers, the tests and `validation/*/README.md`. Reason: a mixed set of samples must say which tool made each row. *Judgement call by the author. Reviewer and user may object.*
- SLURM (`slurm.rs`, `run.rs`): forward `--quantifier`. The required-tool check lists `salmon` instead of `kallisto` when the quantifier is Salmon. Strand-check is unchanged (it runs kallisto or HISAT2 as today; builder confirms it does not depend on the quantifier).
- `mode_word` and the summary `mode` field stay `quant`. The quantifier is a separate field.

### Environment

- `pixi.toml`: add `salmon ==2.8.0` (bioconda). The builder confirms the solve on linux-64. The HPCC module `salmon/1.7.0` is old and is not used.
- Note: the Salmon index is not portable across Salmon versions or architectures. Index per run directory, as kallisto does.

## Tests

- Pure unit tests in core for each builder (paired, single, lanes, three strandedness values, error cases).
- Unit tests for `quant.sf` gene sum (rounding at .5, unknown transcript, empty file, missing column).
- Unit test for the meta parser (good file, missing field).
- Tool-gated test (skips, or fails under `QRNAMAP_REQUIRE_TOOLS=1`, per R5): recover counts on the large synthetic fixture within 10 percent of truth; lanes versus concatenated reads give the same counts; forward and reverse stranded runs.
- Regression test: kallisto gene counts on the large fixture are identical before and after the change.
- CLI test: `--task-id` refuses a marker made with the other quantifier.
- Wasm build check for the core crate still passes.

## Open points for the builder (flags to confirm from `salmon --help` for 2.8.0)

- Exact flag spelling: `--fldMean`, `--fldSD`, `-p`, and any `--ramLimit` for index.
- Whether `quant` needs `--validateMappings` (default in recent versions) or `--minAssignedFrags`.
- Field names in `aux_info/meta_info.json`.
- Whether the output directory must not exist beforehand.

## Validation after the code

Run Salmon on the AMB data (GSE80357) and compare gene counts with the kallisto counts and with the authors' tables, as the GEO comparison did. Report the measured agreement. Do not claim agreement before the numbers exist.

## Decisions after the builder check (author's calls, marked for the user to review)

1. `fraction_assigned` = `fraction_mapped` = `percent_mapped / 100`, the same as kallisto today.
2. Single-end lanes: one space-separated list after `-r`. Paired lanes: all R1 files after `-1`, all R2 files after `-2`. The tool-gated test must show lanes equal concatenated reads.
3. Zero mapped reads: `num_mapped == 0` fails the sample with the existing `NoReadsPseudoaligned` error, as kallisto does. Salmon's `compatible_fragment_ratio` library-type warning is not used.
4. Summary format: the `quantifier` column goes after `mode`. Align rows write `none`. `read_sample_summary` accepts 8 columns only. Old 7-column files are rejected with a clear message.
5. Stale output removal (R24) lists the files of both tools: `abundance.tsv`, `run_info.json`, `quant.sf`, `aux_info/meta_info.json`, and the Salmon per-sample directory. A switch from one tool to the other must leave no old output.
6. `.prepared` marker: add a `quantifier` line, bump the marker header version from 1 to 2, and refuse a marker without the line. The marker stores a value here, not a file stamp. An old prepared directory must be prepared again. The error message says so.
7. `--quantifier salmon` together with align mode is an error. It is not ignored.
8. The Salmon index is a directory. `QuantContext.index` and `require_file` (run.rs) must accept a directory for Salmon.
9. Also in scope, found by the builder: the `check_samples` dry validation branches on the quantifier; kallisto-named error variants and messages (`PairCountMismatch { kallisto }`, `NoReadsPseudoaligned`) take the tool name; `RunFlags::config` sets the quantifier; `slurm_plan` forwards `--quantifier`; `required_tools` takes the quantifier (callers and tests change); help text and README name both tools; tests that hard-code the 7-column header are updated.
10. Minor builder calls: create the sample directory as today; `--quantifier` is a clap `value_enum`; index flags are `-t -i -p` only; add `--validateMappings` only if 2.8.0 needs it; float formatting as kallisto; one shared private sum-and-round helper in `counts.rs`; negative or NaN `NumReads` is a `BadRow` error.
