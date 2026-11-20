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

## Measured

Pipeline recovery on simulated reads (Task 11, 2026-11-19). Tests: `crates/qrnamap-cli/tests/pipeline_quant.rs`, `pipeline_align.rs`.

- Data: 10,000 error-free read pairs from the 6 transcripts of `tests/fixtures/mini` (5 genes, 1,000 to 3,000 pairs per gene). The transcripts are 70 to 120 bases long, shorter than the 200-base fragment, so each fragment is the whole transcript. Reads are 75 bases (70 for the 70-base transcript).
- kallisto 0.52.0 (`-t 1`): every gene count equals the truth (maximum relative error 0). Test tolerance: 10%.
- HISAT2 2.2.2 + featureCounts 2.1.1 (`-p 1`, `-T 1`): every gene count equals the truth (maximum relative error 0); all 10,000 pairs mapped and assigned. Test tolerance: 5%.
- Wrong strandedness: a reverse library run as forward gives featureCounts assigned fraction 0.0 and kallisto 0 pseudoaligned reads (kallisto then exits with status 1). An unstranded library run as forward in kallisto gives a pseudoaligned fraction of 0.496.
- These results come from very short, error-free transcripts. They do not predict accuracy on real data.

Native HISAT2 runs `hisat2-align-s`, not the `hisat2` wrapper: the Perl wrapper passes its arguments through a shell, which split a path with a space and expanded `$HOME` in the tests.

Larger synthetic fixture (`crates/qrnamap-cli/tests/common/large.rs`, test `tests/large.rs`, 2026-11-19): 2 contigs of 20 kb, 12 genes (2 multi-isoform, 1 rRNA, 1 CDS-only), 15 transcripts of 661 to 2,034 bases, 40,000 error-free pairs (read length 100, fragment 200 to 300). The reads come from the generator's own transcripts, not from our builder. Measured maximum relative error per gene (genes with at least 200 pairs): kallisto 0 (unstranded, forward and reverse libraries with the matching flag); HISAT2 + featureCounts 0.0005 unstranded and 0.0018 stranded. The lost pairs are g07 pairs that map to 2 places (`Unassigned_MultiMapping`). A reverse library run as forward: kallisto pseudoaligns no read (exit status 1); featureCounts assigned fraction 0.0.
