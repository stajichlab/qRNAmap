/**
 * Page-v1 measurements: read preparation in a dedicated worker (modes wgz and wplain of measure.mjs, suite quant).
 * Message in: { urls: {1, 2}, pairs, mode }. Message out: { ok, ms, reads: [r1, r2], bytes: [b1, b2] } with the two byte
 * arrays transferred (not copied), or { ok: false, error }.
 *   wgz    qmutfind reads.js keepFirstReads (fflate gunzip, cut at N records, fflate gzip level 3), as qmutfind runs it in its worker.
 *   wplain DecompressionStream('gzip'), cut at N records, plain FASTQ bytes.
 */
import { prepMate, keepFirstPlainDS } from './steps.mjs';

self.onmessage = async ({ data }) => {
  try {
    const started = performance.now();
    const results = [];
    for (const mate of [1, 2]) {
      const response = await fetch(data.urls[mate]);
      if (!response.ok) throw new Error(`fetch mate ${mate}: HTTP ${response.status}`);
      results.push(data.mode === 'wgz'
        ? await prepMate(response.body.getReader(), data.pairs, 'gz', `mate ${mate}`)
        : await keepFirstPlainDS(response.body, data.pairs));
    }
    const bytes = results.map(r => r.bytes);
    self.postMessage({ ok: true, ms: performance.now() - started, reads: results.map(r => r.reads), bytes }, bytes.map(b => b.buffer));
  } catch (error) {
    self.postMessage({ ok: false, error: `${error?.name || 'Error'}: ${error?.message || error}` });
  }
};
