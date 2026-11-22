# Browser page checks

The workflow `.github/workflows/page.yml` does not download or convert reads from NCBI on each run.
It downloads static files from the release `page-v1-test-data-v1` (pre-release) of this repository.

- Default (push, or dispatch with all `counts` <= 1,000,000): the job `reads` downloads the two prefix files, the first 1,000,000 read pairs of each mate of SRR3396381 (83,401,396 and 83,784,618 bytes).
- Full path (the largest value in `counts` > 1,000,000): the job `reads` downloads the two complete mates (712,250,441 and 723,342,896 bytes; 8,405,534 pairs each).
- The SHA-256 of every file is checked against `measure/reads_prefix.tsv`, which is committed in this repository (not against the release). A replaced asset fails the job.
- The prefix files were made from the full files with `zcat SRR3396381_N.fastq.gz | head -n 4000000 | gzip -6 -n`. The first 4,000,000 decompressed lines are identical to those of the full files (SHA-256 compared).
- The full files were made by `measure/prepare_reads.sh` (NCBI SRA, fasterq-dump 3.2.0, gzip -6 -n) in run 44990656616. CI does not run this script. Use it to re-create the data.
- The script `measure/fetch_reads_release.sh <out dir> <prefix|full>` does the download. It also fetches the yeast transcripts from the spike 0 fixture release (SHA-256 from `validation/browser-spike0/manifest.tsv`).
