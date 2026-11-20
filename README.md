# qRNAmap

A Rust command-line tool, with a WASM build, for RNA-seq quantification. Licence: GPL-3.0-only.

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
