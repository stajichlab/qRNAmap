// The preview size N (read pairs per sample) for a read length. N depends on the read length because longer reads need more memory
// (docs/browser/page-v1-measurements.md, "Page quantification path: peak memory"): CI run 45049907640, quant phase at 200,000 pairs,
// whole-browser Pss: 150-base reads (C. neoformans H99) Chromium 1,204 MiB, Firefox 1,663 to 1,759 MiB, WebKit 1,993 to 2,044 MiB,
// above the 1.5 GiB rule; 101-base reads (yeast) up to 1,508 MiB in WebKit.
// 150-base reads, N = 50,000 and 75,000 with the prep worker (CI run 45099586699, 3 runs, 3 browsers): highest Pss of all phases
// WebKit 1,491,448 and 1,467,306 kB (rule: 1,572,864 kB), Firefox 989,096 and 1,094,731 kB, Chromium 645,805 and 666,758 kB.
// N = 100,000 with the prep worker is measured by the next run of page-quant-memory.yml (UNVERIFIED until then). 75,000 is the
// largest value that was measured to pass.
// No imports: the page worker and Node (unit tests) load this file.

/** The largest N the page ever uses (the Task 3 hard limit, measured with 101-base reads only). */
export const HARD_LIMIT = 200000;

/** Rows: reads of at most `maxLength` bases get N = `pairs`. The last row covers every longer read. */
export const PREVIEW_TABLE = [
  { maxLength: 101, pairs: 200000, measured: 'yes: 101-base reads (SRR3396381), Task 3 and run 45049907640' },
  { maxLength: Infinity, pairs: 75000, measured: 'yes: 150-base reads (C. neoformans H99), prep worker, 3 runs in 3 browsers' },
];

/** N for reads of `readLength` bases (the longest first record of the two mates). */
export function previewLimit(readLength) {
  if (!(readLength > 0)) throw new RangeError(`read length must be above 0, not ${readLength}`);
  return PREVIEW_TABLE.find((row) => readLength <= row.maxLength).pairs;
}
