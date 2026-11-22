/**
 * Page-v1 measurements: steps that run the same way in a browser page and in Node.
 * measure.mjs imports this file in the page (/measure/steps.mjs) or in a Node child process.
 * Nothing here reads the network or the disk directly: the caller gives readers (stream.getReader()) and URLs.
 *
 * Reads go to kallisto the way the vendored runtime takes them: whole files in the Emscripten file system (MEMFS).
 * ToolRunner (site/vendor/kallisto/runtime/index.mjs) posts the file bytes to a new worker (structured clone, a copy),
 * and execute.mjs writes them with FS.writeFile (MEMFS keeps its own copy). The caller also keeps its copy until the run ends.
 */
import { keepFirstReads } from './qm/reads.js?v=__BUILD__';
import { VerifyingGunzip } from './qm/gz.js?v=__BUILD__';
import { gunzipSync } from './qm/vendor/fflate.js?v=__BUILD__';

export const KALLISTO_ARGS = {
  // Same order as qrnamap-core toolcall::kallisto_index / kallisto_quant with threads 1, unstranded.
  index: (fasta, idx) => ['index', '-i', idx, fasta],
  quant: (idx, out, m1, m2) => ['quant', '-i', idx, '-o', out, '--plaintext', '-t', '1', m1, m2],
};

export async function sha256(bytes) {
  const digest = await crypto.subtle.digest('SHA-256', bytes);
  return Array.from(new Uint8Array(digest), b => b.toString(16).padStart(2, '0')).join('');
}

export function concat(parts, total) {
  const out = new Uint8Array(total ?? parts.reduce((n, p) => n + p.length, 0));
  let offset = 0;
  for (const p of parts) { out.set(p, offset); offset += p.length; }
  return out;
}

/**
 * Keep the first maxReads FASTQ records (4 lines each) of a gzip stream as plain bytes. Same cut rule as qmutfind reads.js
 * (count newlines, stop at maxReads * 4), but the kept bytes stay plain. A stream that ends before the limit is read to its
 * end and checked (CRC-32, length); an error there is thrown.
 */
export async function keepFirstPlain(reader, maxReads) {
  const target = maxReads * 4;
  const parts = [];
  let total = 0, lines = 0, cut = false, downloaded = 0;
  const feed = bytes => {
    if (cut) return;
    let seen = 0;
    for (let i = 0; i < bytes.length; i++) {
      if (bytes[i] === 10 && ++seen + lines === target) { bytes = bytes.subarray(0, i + 1); cut = true; break; }
    }
    lines += seen;
    parts.push(bytes); // fflate gives a new array per chunk (slc), so no copy is needed
    total += bytes.length;
  };
  const gunzip = new VerifyingGunzip(feed);
  while (!cut) {
    const { done, value } = await reader.read();
    if (done) { gunzip.push(new Uint8Array(0), true); break; }
    downloaded += value.length;
    try { gunzip.push(value, false); } catch (error) { if (!cut) throw error; }
  }
  if (cut) await reader.cancel();
  return { bytes: concat(parts, total), reads: Math.floor(lines / 4), cut, downloaded };
}

/** Same cut rule, with the built-in DecompressionStream('gzip') instead of fflate. Gives plain bytes. */
export async function keepFirstPlainDS(stream, maxReads) {
  const reader = stream.pipeThrough(new DecompressionStream('gzip')).getReader();
  const target = maxReads * 4;
  const parts = [];
  let total = 0, lines = 0, cut = false;
  while (!cut) {
    const { done, value } = await reader.read();
    if (done) break;
    let bytes = value, seen = 0;
    for (let i = 0; i < bytes.length; i++) {
      if (bytes[i] === 10 && ++seen + lines === target) { bytes = bytes.subarray(0, i + 1); cut = true; break; }
    }
    lines += seen;
    parts.push(bytes);
    total += bytes.length;
  }
  if (cut) await reader.cancel().catch(() => {});
  return { bytes: concat(parts, total), reads: Math.floor(lines / 4), cut };
}

/** Prepare one mate: mode 'gz' uses qmutfind keepFirstReads (fflate gzip level 3 of the kept reads); 'plain' keeps plain bytes. */
export async function prepMate(reader, maxReads, mode, label) {
  if (mode === 'gz') {
    const r = await keepFirstReads(reader, { maxReads, maxKept: 2 ** 40, maxDownload: 2 ** 40, label, progressEvery: 0 });
    return { bytes: r.bytes, reads: r.reads, cut: r.cut };
  }
  return keepFirstPlain(reader, maxReads);
}

/** Read names (first word of each header line, without a /1 or /2 end) of a plain FASTQ. */
export function recordNames(plain) {
  const names = [];
  const text = new TextDecoder();
  let start = 0, line = 0;
  for (let i = 0; i < plain.length; i++) {
    if (plain[i] !== 10) continue;
    if (line % 4 === 0) {
      let end = start;
      while (end < i && plain[end] !== 32 && plain[end] !== 9) end++;
      names.push(text.decode(plain.subarray(start + 1, end)).replace(/\/[12]$/, ''));
    }
    line++;
    start = i + 1;
  }
  return names;
}

/** Compare the record names of two mates, record by record. */
export function pairCheck(plain1, plain2) {
  const a = recordNames(plain1), b = recordNames(plain2);
  let firstDiff = -1;
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) if (a[i] !== b[i]) { firstDiff = i; break; }
  return { records1: a.length, records2: b.length, names_equal: firstDiff === -1 && a.length === b.length, first_diff: firstDiff };
}

export function plainOf(bytes) {
  return bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes;
}

/**
 * One kallisto command through the vendored ToolRunner. files: [{path, data}]. Gives exit code, stderr text, output files,
 * heap record (from heap_wrap.mjs), time in ms, or an error string.
 */
export async function runKallisto({ runtimeUrl, kallistoUrl, wrapUrl, args, files, outputs, timeoutMs }) {
  const { ToolRunner } = await import(runtimeUrl);
  let broadcast = null;
  const channel = new BroadcastChannel('measure-heap');
  channel.onmessage = event => { broadcast = event.data; };
  const moduleUrl = new URL(`${wrapUrl}?tool=${encodeURIComponent(kallistoUrl)}`);
  const runner = new ToolRunner({ moduleUrl, waitForExit: false });
  const started = performance.now();
  const result = {};
  try {
    const run = await runner.run(args, { files, outputs: [...outputs, 'measure_heap.json'], timeoutMs });
    result.ms = performance.now() - started;
    result.exitCode = run.exitCode;
    result.missing = run.missingOutputs.filter(p => !p.endsWith('measure_heap.json'));
    result.files = {};
    for (const file of run.files) {
      const name = file.path.replace(/^\/work\//, '');
      if (name === 'measure_heap.json') { result.heap = JSON.parse(new TextDecoder().decode(file.data)); continue; }
      result.files[name] = file.data;
    }
    result.stderr = new TextDecoder().decode(run.stderr);
  } catch (error) {
    result.ms = performance.now() - started;
    result.error = `${error?.name || 'Error'}: ${error?.message || error}`.slice(0, 1200);
  }
  await new Promise(done => setTimeout(done, 200));
  channel.close();
  result.heap = result.heap || broadcast;
  return result;
}

/** Facts from kallisto run_info.json (n_processed, n_pseudoaligned, p_pseudoaligned). */
export function runInfo(bytes) {
  try {
    const info = JSON.parse(new TextDecoder().decode(bytes));
    return { n_processed: info.n_processed, n_pseudoaligned: info.n_pseudoaligned, p_pseudoaligned: info.p_pseudoaligned };
  } catch { return {}; }
}
