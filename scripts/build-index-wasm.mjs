// Build a kallisto index with the kallisto.wasm build the page uses (site/vendor/kallisto), in Node.
//   node scripts/build-index-wasm.mjs TRANSCRIPTS.fa[.gz] OUT.idx
// The index file format differs between kallisto builds of the same version: an index from the Linux (conda) kallisto
// 0.52.0 fails in kallisto.wasm 0.52.0 and in the macOS Homebrew build with "Corrupted index; unitig not found"
// (bench 2026-11-26_index-compat). Indexes for the page must therefore come from kallisto.wasm itself. Same argument
// list as the page (kallisto_args kind "index": k = 31). Prints the build time and the index size.
import { readFileSync, writeFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { dirname, resolve } from 'node:path';
import { execute } from '../site/vendor/kallisto/runtime/execute.mjs';

const [fasta, out] = process.argv.slice(2);
if (!fasta || !out) {
  console.error('usage: node scripts/build-index-wasm.mjs TRANSCRIPTS.fa[.gz] OUT.idx');
  process.exit(2);
}
const here = dirname(fileURLToPath(import.meta.url));
const moduleUrl = pathToFileURL(resolve(here, '../site/vendor/kallisto/kallisto.mjs')).href;
let data = readFileSync(fasta);
if (data[0] === 0x1f && data[1] === 0x8b) data = gunzipSync(data);
const started = performance.now();
const result = await execute({
  moduleUrl,
  args: ['index', '-i', 'transcripts.idx', '-k', '31', 'transcripts.fa'],
  files: [{ path: '/work/transcripts.fa', data }],
  outputs: ['/work/transcripts.idx'],
});
data = null;
if (result.exitCode !== 0 || result.missingOutputs.length) {
  console.error(new TextDecoder().decode(result.stderr).split('\n').slice(-5).join('\n'));
  console.error(`kallisto.wasm index failed (exit ${result.exitCode})`);
  process.exit(1);
}
writeFileSync(out, result.files[0].data);
console.log(`built ${out}: ${result.files[0].data.length} bytes in ${((performance.now() - started) / 1000).toFixed(1)} s`);
