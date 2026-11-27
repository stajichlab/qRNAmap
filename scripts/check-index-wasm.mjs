// Check that a kallisto index works in the kallisto.wasm build the page uses (site/vendor/kallisto), in Node.
//   node scripts/check-index-wasm.mjs INDEX[.gz] R1.fastq[.gz] [R2.fastq[.gz]]
// Runs `kallisto quant` (threads 1, as the page) with the index and the reads, prints one line
// `ok <processed> <pseudoaligned>` or `FAIL <reason>` and exits 0 or 1. A gzip index (as in the release) is gunzipped first.
// Used by .github/workflows/build-indexes.yml on every built index before release, so an index the page cannot read
// is never published. The reads can be any small FASTQ (validation/browser-page/fixtures/s1_R*.fastq.gz).
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, resolve } from 'node:path';
import { execute } from '../site/vendor/kallisto/runtime/execute.mjs';

const [indexPath, r1, r2] = process.argv.slice(2);
if (!indexPath || !r1) {
  console.error('usage: node scripts/check-index-wasm.mjs INDEX[.gz] R1.fastq[.gz] [R2.fastq[.gz]]');
  process.exit(2);
}
const here = dirname(fileURLToPath(import.meta.url));
const moduleUrl = pathToFileURL(resolve(here, '../site/vendor/kallisto/kallisto.mjs')).href;
let index = readFileSync(indexPath);
if (index[0] === 0x1f && index[1] === 0x8b) index = gunzipSync(index);
const reads = [r1, r2].filter(Boolean);
const names = reads.map((p, i) => `r${i + 1}.fastq${p.endsWith('.gz') ? '.gz' : ''}`);
const args = ['quant', '-i', 'i.idx', '-o', 'kq', '--plaintext', '-t', '1',
  ...(reads.length === 1 ? ['--single', '-l', '200', '-s', '20'] : []), ...names];
const result = await execute({
  moduleUrl, args,
  files: [{ path: '/work/i.idx', data: index }, ...reads.map((p, i) => ({ path: `/work/${names[i]}`, data: readFileSync(p) }))],
  outputs: ['/work/kq/run_info.json'],
});
const stderr = new TextDecoder().decode(result.stderr);
if (result.exitCode !== 0 || result.missingOutputs.length) {
  const why = stderr.split('\n').filter((l) => /Error|error/.test(l)).join(' ') || `exit ${result.exitCode}`;
  console.log(`FAIL ${why}`);
  process.exit(1);
}
const info = JSON.parse(new TextDecoder().decode(result.files[0].data));
console.log(`ok ${info.n_processed} ${info.n_pseudoaligned}`);
