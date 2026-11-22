#!/usr/bin/env node
/**
 * Page-v1 measurements (Task 3 of docs/superpowers/plans/2026-11-22-browser-page-v1.md).
 * Method after validation/browser-spike0/index_size.mjs: Node or a headless Playwright browser, one fresh process or browser
 * launch per run, whole-browser Pss from /proc, Node maxRSS, WASM heap from heap_wrap.mjs, a local static server.
 * Tools: the vendored kallisto (site/vendor/kallisto) and its ToolRunner. Reads and genomes: files in --data.
 *
 * Usage (from the repository root):
 *   node validation/browser-page/measure/measure.mjs --suite quant|truncated|genomes --engine node|chromium|firefox|webkit
 *        --data DIR --out DIR [--counts 100000,200000] [--modes gz,plain] [--reps 3] [--pkg DIR] [--genomes FILE]
 *
 * Suites:
 *   quant      per count and mode: kallisto index of yeast.transcripts.fa, then the first COUNT pairs of
 *              SRR3396381_1/_2.fastq.gz (prepared in the page: mode gz = qmutfind reads.js keepFirstReads, which packs the
 *              kept reads as gzip level 3; mode plain = the same cut, kept as plain FASTQ), then kallisto quant.
 *   truncated  browsers only: Range requests to ENA, gzip decoders on a cut stream, kallisto on cut files, the N-pair method.
 *   genomes    per genome (--genomes TSV: accession, fna file, gff file, fna URL, gff URL): gunzip (DecompressionStream and
 *              fflate), build_transcriptome (qrnamap-wasm), kallisto index.
 * Probes (see the spike README for their limits): tree_pss (sum of Pss over all processes below this Node process, every
 * 250 ms), node_rss (child RSS every 250 ms), node_maxrss_kb (kernel peak of the child), wasm heap (HEAPU8.length after callMain).
 * Guards: a run is stopped when the tree Pss or the child RSS passes --guard-gib (default 10 GiB of the 16 GB runner); the
 * outcome is then mem_guard. A count is the last one tried when no run of it finished, or when a run passed --stop-gib
 * (default 6 GiB); larger counts get outcome skipped. A failed run is a result; the next run still starts.
 */
import { readFile, writeFile, mkdir, stat } from 'node:fs/promises';
import { readFileSync, readdirSync, createReadStream } from 'node:fs';
import { createHash } from 'node:crypto';
import { createServer } from 'node:http';
import { spawn } from 'node:child_process';
import { resolve, join, sep, dirname } from 'node:path';
import { pathToFileURL, fileURLToPath } from 'node:url';
import { Readable } from 'node:stream';
import { createRequire } from 'node:module';

const options = { suite: 'quant', engine: 'node', data: '', out: '', root: process.cwd(), counts: '100000,200000,500000,1000000',
  modes: 'gz,plain', reps: '3', pkg: '', genomes: '', 'guard-gib': '10', 'stop-gib': '6', 'timeout-min': '30', 'budget-min': '300',
  child: '', count: '', mode: '', rep: '', genome: '', 'n-method': '50000', 'range-mib': '8' };
for (let i = 2; i < process.argv.length; i += 2) {
  const key = process.argv[i].replace(/^--/, '');
  if (!(key in options) || process.argv[i + 1] === undefined) throw new Error(`unknown or empty option ${process.argv[i]}`);
  options[key] = process.argv[i + 1];
}
const ENGINES = ['node', 'chromium', 'firefox', 'webkit'];
if (!ENGINES.includes(options.engine)) throw new Error(`engine must be one of ${ENGINES}`);
if (!options.data || !options.out) throw new Error('--data and --out are required');
const root = resolve(options.root), dataDir = resolve(options.data), outDir = resolve(options.out);
const here = dirname(fileURLToPath(import.meta.url));
const vendorDir = join(root, 'site', 'vendor', 'kallisto');
const pkgDir = options.pkg ? resolve(options.pkg) : '';
const SAMPLE_MS = 250;
const TIMEOUT_MS = Number(options['timeout-min']) * 60000;
const GUARD_KB = Number(options['guard-gib']) * 1024 * 1024;
const STOP_KB = Number(options['stop-gib']) * 1024 * 1024;
const BUDGET_MS = Number(options['budget-min']) * 60000;
const RUN = 'SRR3396381';
const ENA = { 1: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/001/SRR3396381/SRR3396381_1.fastq.gz',
  2: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/001/SRR3396381/SRR3396381_2.fastq.gz' };
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const log = message => console.error(`[${new Date().toISOString()}] ${options.engine} ${message}`);
const kb = bytes => Math.round(bytes / 1024);
const s2 = ms => (ms / 1000).toFixed(2);
const clean = value => String(value ?? '').replace(/[\t\r\n]+/g, ' ');

function errorOutcome(text) {
  if (/mem_guard/.test(text)) return 'mem_guard';
  if (/\bOOM\b|out of memory|Cannot enlarge memory|bad_alloc|Array buffer allocation failed|allocation fail|RangeError: .*(memory|length|size)|Invalid typed array length|Invalid array buffer length/i.test(text)) return 'oom';
  if (/crash|Target (page, context or browser )?closed|has been closed|disconnected|exited without a result/i.test(text)) return 'crash';
  if (/exceeded \d+ ms|TimeoutError/.test(text)) return 'timeout';
  return 'error';
}

// ---------------------------------------------------------------- process tree sampling (browsers)
function readdirSafe(path) { try { return readdirSync(path); } catch { return []; } }
function procTree(rootPid = process.pid) {
  const parent = new Map();
  for (const name of readdirSafe('/proc')) {
    if (!/^\d+$/.test(name)) continue;
    try {
      const st = readFileSync(`/proc/${name}/stat`, 'utf8');
      parent.set(Number(name), Number(st.slice(st.lastIndexOf(')') + 2).split(' ')[1]));
    } catch {}
  }
  const out = [];
  for (const pid of parent.keys()) {
    let p = parent.get(pid), depth = 0;
    while (p && p !== rootPid && depth++ < 64) p = parent.get(p);
    if (p === rootPid) out.push(pid);
  }
  return out;
}
function sampleTree() {
  let pss = 0, rss = 0, maxRss = 0, procs = 0, pssMissing = 0;
  for (const pid of procTree()) {
    let p = null, r = null;
    try {
      const text = readFileSync(`/proc/${pid}/smaps_rollup`, 'utf8');
      p = Number(/^Pss:\s+(\d+) kB/m.exec(text)?.[1]);
      r = Number(/^Rss:\s+(\d+) kB/m.exec(text)?.[1]);
    } catch {}
    if (r === null || Number.isNaN(r)) {
      try { r = Number(/^VmRSS:\s+(\d+) kB/m.exec(readFileSync(`/proc/${pid}/status`, 'utf8'))?.[1]); } catch { continue; }
    }
    if (Number.isNaN(r)) continue;
    procs++;
    if (p === null || Number.isNaN(p)) { pssMissing++; p = r; }
    pss += p; rss += r; maxRss = Math.max(maxRss, r);
  }
  return { pss, rss, maxRss, procs, pssMissing };
}
class Sampler {
  constructor(onGuard) { this.onGuard = onGuard; this.reset(); this.timer = setInterval(() => this.sample(), SAMPLE_MS); }
  reset() { this.peak = { pss: 0, rss: 0, maxRss: 0, procs: 0, pssMissing: 0, samples: 0 }; this.sample(); }
  sample() {
    const s = sampleTree();
    const p = this.peak;
    p.samples++; p.pss = Math.max(p.pss, s.pss); p.rss = Math.max(p.rss, s.rss); p.maxRss = Math.max(p.maxRss, s.maxRss);
    p.procs = Math.max(p.procs, s.procs); p.pssMissing = Math.max(p.pssMissing, s.pssMissing);
    if (s.pss > GUARD_KB && !this.tripped) { this.tripped = s.pss; this.onGuard(s.pss); }
    return s;
  }
  stop() { clearInterval(this.timer); }
}
function killTree() { for (const pid of procTree()) { try { process.kill(pid, 'SIGKILL'); } catch {} } }

// ---------------------------------------------------------------- static server
let server;
function contentType(path) {
  if (path.endsWith('.wasm')) return 'application/wasm';
  if (/\.(m?js)$/.test(path)) return 'text/javascript';
  if (path.endsWith('.html')) return 'text/html';
  return 'application/octet-stream';
}
async function sendFile(response, base, rel) {
  const path = resolve(base, '.' + rel);
  if (path !== base && !path.startsWith(base + sep)) { response.writeHead(403).end(); return; }
  const info = await stat(path);
  response.writeHead(200, { 'Content-Type': contentType(path), 'Content-Length': info.size });
  createReadStream(path).pipe(response);
}
function startServer() {
  server = createServer(async (request, response) => {
    try {
      const name = decodeURIComponent(new URL(request.url, 'http://localhost').pathname); // the query (?v=...) is ignored
      if (name === '/') { response.setHeader('Content-Type', 'text/html'); response.end('<!doctype html><meta charset="utf-8"><title>page-v1 measure</title>'); return; }
      if (name.startsWith('/kallisto/')) return await sendFile(response, vendorDir, name.slice('/kallisto'.length));
      if (name.startsWith('/measure/')) return await sendFile(response, here, name.slice('/measure'.length));
      if (name.startsWith('/pkg/') && pkgDir) return await sendFile(response, pkgDir, name.slice('/pkg'.length));
      if (name.startsWith('/data/')) return await sendFile(response, dataDir, name.slice('/data'.length));
      response.writeHead(404).end();
    } catch (error) { log(`server ${request.url}: ${error.message}`); if (!response.headersSent) response.writeHead(404); response.end(); }
  });
  return new Promise(done => server.listen(0, '127.0.0.1', done));
}

// ---------------------------------------------------------------- browser session
async function openBrowser(playwright) {
  const browser = await playwright[options.engine].launch({ headless: true });
  const page = await browser.newPage();
  page.on('console', m => { if (m.type() !== 'debug') log(`browser ${m.type()}: ${m.text().slice(0, 300)}`); });
  page.on('pageerror', e => log(`browser pageerror: ${e.message}`));
  let fail;
  const failed = new Promise((_, reject) => { fail = reject; });
  failed.catch(() => {});
  page.on('crash', () => fail(new Error('page crash event')));
  browser.on('disconnected', () => fail(new Error('browser disconnected')));
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  const session = { browser, page, failed, fail, version: browser.version() };
  session.sampler = new Sampler(pss => { log(`mem_guard: tree Pss ${pss} kB`); fail(new Error(`mem_guard: tree Pss ${pss} kB passed ${GUARD_KB} kB`)); killTree(); });
  session.eval = async (fn, arg) => {
    let guard;
    try {
      return await Promise.race([page.evaluate(fn, arg), failed, new Promise((_, reject) => {
        guard = setTimeout(() => reject(new Error(`browser run exceeded ${TIMEOUT_MS} ms`)), TIMEOUT_MS);
      })]);
    } finally { clearTimeout(guard); }
  };
  return session;
}
async function closeBrowser(session) {
  if (!session) return;
  session.sampler?.stop();
  await Promise.race([session.browser.close().catch(() => {}), new Promise(done => setTimeout(done, 15000))]);
  killTree();
}
const peakCols = (prefix, p) => ({ [`${prefix}_tree_pss_peak_kb`]: p.pss, [`${prefix}_tree_rss_sum_peak_kb`]: p.rss, [`${prefix}_samples`]: p.samples });

// URLs of the tools in the page and in Node
const pageUrls = () => ({ runtimeUrl: '/kallisto/runtime/index.mjs', kallistoUrl: new URL('/kallisto/kallisto.mjs', `http://127.0.0.1:${server.address().port}`).href,
  wrapUrl: new URL('/measure/heap_wrap.mjs', `http://127.0.0.1:${server.address().port}`).href });
const nodeUrls = () => ({ runtimeUrl: pathToFileURL(join(vendorDir, 'runtime', 'index.mjs')).href, kallistoUrl: pathToFileURL(join(vendorDir, 'kallisto.mjs')).href,
  wrapUrl: pathToFileURL(join(here, 'heap_wrap.mjs')).href });

// ================================================================= suite quant
const QUANT_COLS = ['engine', 'suite', 'mode', 'pairs', 'rep', 'outcome', 'phase_failed', 'index_s', 'prep_s', 'quant_s', 'total_s',
  'reads1', 'reads2', 'input1_bytes', 'input2_bytes', 'n_processed', 'n_pseudoaligned', 'p_pseudoaligned',
  'index_heap_final_bytes', 'quant_heap_initial_bytes', 'quant_heap_final_bytes', 'quant_callmain_s',
  'index_tree_pss_peak_kb', 'prep_tree_pss_peak_kb', 'quant_tree_pss_baseline_kb', 'quant_tree_pss_peak_kb', 'quant_tree_rss_sum_peak_kb',
  'tree_pss_load_kb', 'node_maxrss_kb', 'node_rss_index_peak_kb', 'node_rss_prep_peak_kb', 'node_rss_quant_peak_kb',
  'abundance_sha256', 'index_bytes', 'browser_version', 'detail'];

// In-page (or in-Node) phases. These functions are serialised by Playwright, so they import steps.mjs themselves.
async function phaseIndex(cfg) {
  const s = await import(cfg.stepsUrl);
  let tx = globalThis.measureTx;
  if (!tx) { const r = await fetch(cfg.txUrl); if (!r.ok) throw new Error(`fetch ${cfg.txUrl}: HTTP ${r.status}`); tx = new Uint8Array(await r.arrayBuffer()); }
  const r = await s.runKallisto({ ...cfg.urls, args: s.KALLISTO_ARGS.index('yeast.transcripts.fa', 'yeast.idx'),
    files: [{ path: 'yeast.transcripts.fa', data: tx }], outputs: ['yeast.idx'], timeoutMs: cfg.timeoutMs });
  globalThis.measureTx = null;
  globalThis.measureIdx = r.files?.['yeast.idx'] || null;
  return { ms: r.ms, exitCode: r.exitCode, error: r.error, missing: r.missing, heap: r.heap, indexBytes: globalThis.measureIdx?.length ?? 0,
    stderrTail: (r.stderr || '').slice(-300) };
}
async function phasePrep(cfg) {
  const s = await import(cfg.stepsUrl);
  const started = performance.now();
  const out = { };
  globalThis.measureReads = {};
  for (const mate of [1, 2]) {
    const reader = cfg.open ? await cfg.open(mate) : (await (async () => {
      const r = await fetch(cfg.readUrls[mate]); if (!r.ok) throw new Error(`fetch mate ${mate}: HTTP ${r.status}`); return r.body.getReader();
    })());
    const m = await s.prepMate(reader, cfg.pairs, cfg.mode, `mate ${mate}`);
    globalThis.measureReads[mate] = m.bytes;
    out[`reads${mate}`] = m.reads; out[`bytes${mate}`] = m.bytes.length;
  }
  out.ms = performance.now() - started;
  return out;
}
async function phaseQuant(cfg) {
  const s = await import(cfg.stepsUrl);
  const ext = cfg.mode === 'gz' ? '.fastq.gz' : '.fastq';
  const files = [{ path: 'yeast.idx', data: globalThis.measureIdx }, { path: `r_1${ext}`, data: globalThis.measureReads[1] }, { path: `r_2${ext}`, data: globalThis.measureReads[2] }];
  const r = await s.runKallisto({ ...cfg.urls, args: s.KALLISTO_ARGS.quant('yeast.idx', 'kq', `r_1${ext}`, `r_2${ext}`), files,
    outputs: ['kq/abundance.tsv', 'kq/run_info.json'], timeoutMs: cfg.timeoutMs });
  const out = { ms: r.ms, exitCode: r.exitCode, error: r.error, missing: r.missing, heap: r.heap, stderrTail: (r.stderr || '').slice(-400) };
  if (r.files?.['kq/abundance.tsv']) out.abundance = await s.sha256(r.files['kq/abundance.tsv']);
  if (r.files?.['kq/run_info.json']) out.info = s.runInfo(r.files['kq/run_info.json']);
  globalThis.measureReads = null;
  return out;
}

function fillPhase(row, name, r) {
  if (r.error || r.exitCode !== 0 || (r.missing && r.missing.length)) {
    row.phase_failed = name;
    row.detail = r.error || `exit ${r.exitCode}; missing ${(r.missing || []).join(',') || 'none'}; stderr tail: ${r.stderrTail}`;
    row.outcome = r.error ? errorOutcome(r.error) : 'fail';
    return false;
  }
  return true;
}

async function quantBrowserRun(playwright, pairs, mode, rep) {
  const row = { engine: options.engine, suite: 'quant', mode, pairs, rep };
  let session;
  try {
    session = await openBrowser(playwright);
    row.browser_version = session.version;
    const port = server.address().port;
    const cfg = { stepsUrl: `http://127.0.0.1:${port}/measure/steps.mjs`, urls: pageUrls(), timeoutMs: TIMEOUT_MS, pairs, mode,
      txUrl: '/data/yeast.transcripts.fa', readUrls: { 1: `/data/${RUN}_1.fastq.gz`, 2: `/data/${RUN}_2.fastq.gz` } };
    await new Promise(done => setTimeout(done, 1000));
    row.tree_pss_load_kb = sampleTree().pss;
    session.sampler.reset();
    const ix = await session.eval(phaseIndex, cfg);
    Object.assign(row, { index_s: s2(ix.ms), index_heap_final_bytes: ix.heap?.finalBytes, index_bytes: ix.indexBytes, index_tree_pss_peak_kb: session.sampler.peak.pss });
    if (!fillPhase(row, 'index', ix)) return row;
    session.sampler.reset();
    const pr = await session.eval(phasePrep, cfg);
    Object.assign(row, { prep_s: s2(pr.ms), reads1: pr.reads1, reads2: pr.reads2, input1_bytes: pr.bytes1, input2_bytes: pr.bytes2, prep_tree_pss_peak_kb: session.sampler.peak.pss });
    await new Promise(done => setTimeout(done, 500));
    session.sampler.reset();
    row.quant_tree_pss_baseline_kb = session.sampler.peak.pss;
    const q = await session.eval(phaseQuant, cfg);
    Object.assign(row, { quant_s: s2(q.ms), quant_heap_initial_bytes: q.heap?.initialBytes, quant_heap_final_bytes: q.heap?.finalBytes,
      quant_callmain_s: q.heap?.callMainMs != null ? s2(q.heap.callMainMs) : '', quant_tree_pss_peak_kb: session.sampler.peak.pss,
      quant_tree_rss_sum_peak_kb: session.sampler.peak.rss, abundance_sha256: q.abundance, ...(q.info || {}) });
    if (!fillPhase(row, 'quant', q)) return row;
    row.total_s = s2(ix.ms + pr.ms + q.ms);
    row.outcome = Number(row.n_processed) === pairs && row.reads1 === pairs && row.reads2 === pairs ? 'finished' : 'fail';
    if (row.outcome === 'fail') row.detail = `n_processed ${row.n_processed} or record counts ${row.reads1}/${row.reads2} differ from ${pairs}`;
  } catch (error) {
    row.detail = String(error?.message || error).slice(0, 1200);
    row.outcome = errorOutcome(row.detail);
    if (session?.sampler) row.quant_tree_pss_peak_kb = row.quant_tree_pss_peak_kb ?? session.sampler.peak.pss;
    row.phase_failed = row.phase_failed || (row.prep_s ? 'quant' : row.index_s ? 'prep' : 'index');
  } finally {
    await closeBrowser(session);
  }
  return row;
}

// Node child: one quant run in this process.
async function quantChild() {
  const pairs = Number(options.count), mode = options.mode;
  const row = { engine: 'node', suite: 'quant', mode, pairs, rep: options.rep };
  let rssPeak = 0, samples = 0;
  const timer = setInterval(() => { samples++; rssPeak = Math.max(rssPeak, process.memoryUsage().rss); }, SAMPLE_MS);
  const phase = () => { const p = Math.max(rssPeak, process.memoryUsage().rss); rssPeak = 0; return kb(p); };
  const cfg = { stepsUrl: pathToFileURL(join(here, 'steps.mjs')).href, urls: nodeUrls(), timeoutMs: TIMEOUT_MS, pairs, mode,
    open: async mate => Readable.toWeb(createReadStream(join(dataDir, `${RUN}_${mate}.fastq.gz`))).getReader() };
  try {
    globalThis.measureTx = new Uint8Array(await readFile(join(dataDir, 'yeast.transcripts.fa')));
    phase();
    const ix = await phaseIndex(cfg);
    Object.assign(row, { index_s: s2(ix.ms), index_heap_final_bytes: ix.heap?.finalBytes, index_bytes: ix.indexBytes, node_rss_index_peak_kb: phase() });
    if (fillPhase(row, 'index', ix)) {
      const pr = await phasePrep(cfg);
      Object.assign(row, { prep_s: s2(pr.ms), reads1: pr.reads1, reads2: pr.reads2, input1_bytes: pr.bytes1, input2_bytes: pr.bytes2, node_rss_prep_peak_kb: phase() });
      const q = await phaseQuant(cfg);
      Object.assign(row, { quant_s: s2(q.ms), quant_heap_initial_bytes: q.heap?.initialBytes, quant_heap_final_bytes: q.heap?.finalBytes,
        quant_callmain_s: q.heap?.callMainMs != null ? s2(q.heap.callMainMs) : '', node_rss_quant_peak_kb: phase(), abundance_sha256: q.abundance, ...(q.info || {}) });
      if (fillPhase(row, 'quant', q)) {
        row.total_s = s2(ix.ms + pr.ms + q.ms);
        row.outcome = Number(row.n_processed) === pairs && row.reads1 === pairs && row.reads2 === pairs ? 'finished' : 'fail';
        if (row.outcome === 'fail') row.detail = `n_processed ${row.n_processed} or record counts ${row.reads1}/${row.reads2} differ from ${pairs}`;
      }
    }
  } catch (error) {
    row.detail = String(error?.stack || error).slice(0, 1200);
    row.outcome = errorOutcome(row.detail);
  }
  clearInterval(timer);
  row.node_maxrss_kb = process.resourceUsage().maxRSS;
  row.samples = samples;
  process.stdout.write('MEASURE_RESULT ' + JSON.stringify(row) + '\n');
}

function runNodeChild(extra, fallback) {
  return new Promise(done => {
    const args = [fileURLToPath(import.meta.url), '--suite', options.suite, '--engine', 'node', '--data', dataDir, '--out', outDir, '--root', root,
      '--timeout-min', options['timeout-min'], ...(pkgDir ? ['--pkg', pkgDir] : []), ...extra];
    const child = spawn(process.execPath, args, { stdio: ['ignore', 'pipe', 'pipe'] });
    let stdout = '', stderr = '', guarded = 0;
    const watch = setInterval(() => {
      try {
        const rss = Number(/^VmRSS:\s+(\d+) kB/m.exec(readFileSync(`/proc/${child.pid}/status`, 'utf8'))?.[1]);
        if (rss > GUARD_KB && !guarded) { guarded = rss; log(`mem_guard: child RSS ${rss} kB`); child.kill('SIGKILL'); }
      } catch {}
    }, SAMPLE_MS);
    const timer = setTimeout(() => child.kill('SIGKILL'), TIMEOUT_MS + 120000);
    child.stdout.on('data', c => { stdout += c; });
    child.stderr.on('data', c => { stderr += c; process.stderr.write(c); });
    child.on('close', (code, signal) => {
      clearInterval(watch); clearTimeout(timer);
      const line = stdout.split('\n').find(l => l.startsWith('MEASURE_RESULT '));
      if (line) { done(JSON.parse(line.slice(15))); return; }
      const detail = guarded ? `mem_guard: child RSS ${guarded} kB passed ${GUARD_KB} kB` : `child exit ${code} signal ${signal}; stderr tail: ${stderr.slice(-800)}`;
      done({ ...fallback, outcome: guarded ? 'mem_guard' : (signal || code === 134 || code === 137 ? 'crash' : errorOutcome(detail)), detail });
    });
  });
}

function toTsv(cols, row) { return cols.map(k => clean(row[k])).join('\t'); }

async function suiteQuant(playwright, meta) {
  const counts = options.counts.split(',').map(Number), modes = options.modes.split(','), reps = Number(options.reps);
  const started = Date.now();
  for (const mode of modes) {
    const tsv = join(outDir, `quant.${options.engine}.${mode}.tsv`);
    const lines = [QUANT_COLS.join('\t')];
    let stopped = '';
    for (const pairs of counts) {
      let finished = 0, maxPeak = 0;
      for (let rep = 1; rep <= reps; rep++) {
        let row;
        if (stopped) row = { engine: options.engine, suite: 'quant', mode, pairs, rep, outcome: 'skipped', detail: stopped };
        else if (Date.now() - started > BUDGET_MS) row = { engine: options.engine, suite: 'quant', mode, pairs, rep, outcome: 'skipped', detail: `job time budget of ${options['budget-min']} min used` };
        else {
          log(`quant ${mode} ${pairs} pairs rep ${rep} starting`);
          row = options.engine === 'node'
            ? await runNodeChild(['--child', 'quant', '--count', String(pairs), '--mode', mode, '--rep', String(rep)], { engine: 'node', suite: 'quant', mode, pairs, rep })
            : await quantBrowserRun(playwright, pairs, mode, rep);
        }
        if (row.outcome === 'finished') finished++;
        maxPeak = Math.max(maxPeak, Number(row.quant_tree_pss_peak_kb || 0), Number(row.node_maxrss_kb || 0));
        lines.push(toTsv(QUANT_COLS, row));
        await writeFile(tsv, lines.join('\n') + '\n');
        log(`quant ${mode} ${pairs} rep ${rep}: ${row.outcome} total_s=${row.total_s ?? ''} quant_s=${row.quant_s ?? ''} peak=${row.quant_tree_pss_peak_kb ?? row.node_maxrss_kb ?? ''} ${row.detail ?? ''}`);
      }
      if (!stopped && finished === 0) stopped = `no run finished at ${pairs} pairs`;
      if (!stopped && maxPeak > STOP_KB) stopped = `a run at ${pairs} pairs passed ${options['stop-gib']} GiB (${maxPeak} kB)`;
    }
  }
}

// ================================================================= suite truncated (browsers)
async function truncatedInPage(cfg) {
  const s = await import(cfg.stepsUrl);
  const rows = [];
  const add = (test, result, detail) => rows.push({ test, result, detail });
  const errText = e => `${e?.name || 'Error'}: ${e?.message || e}`.slice(0, 500);
  const range = async (url, bytes) => {
    const t = performance.now();
    const r = await fetch(url, { headers: { Range: `bytes=0-${bytes - 1}` }, mode: 'cors', credentials: 'omit' });
    const body = new Uint8Array(await r.arrayBuffer());
    return { status: r.status, body, ms: performance.now() - t, contentRange: r.headers.get('Content-Range') };
  };
  const R = cfg.rangeMiB * 1048576, R2 = Math.round(cfg.rangeMiB * 0.75) * 1048576;
  const cut = {};
  for (const mate of [1, 2]) {
    try {
      cut[mate] = await range(cfg.ena[mate], R);
      add(`range_fetch_mate${mate}`, cut[mate].status === 206 && cut[mate].body.length === R ? 'pass' : 'fail',
        `status=${cut[mate].status}; bytes=${cut[mate].body.length}; ms=${cut[mate].ms.toFixed(0)}; kB_per_s=${(cut[mate].body.length / cut[mate].ms).toFixed(1)}; magic=${cut[mate].body[0]?.toString(16)}${cut[mate].body[1]?.toString(16)}; content_range_readable=${cut[mate].contentRange ?? 'null'}`);
    } catch (e) { add(`range_fetch_mate${mate}`, 'fail', errText(e)); }
  }
  if (!cut[1]?.body?.length || !cut[2]?.body?.length) return rows;
  // Decoders on a cut stream.
  for (const mate of [1, 2]) {
    let plain = 0, lines = 0, error = 'none';
    try {
      const reader = new Blob([cut[mate].body]).stream().pipeThrough(new DecompressionStream('gzip')).getReader();
      for (;;) { const { done, value } = await reader.read(); if (done) break; plain += value.length; for (const b of value) if (b === 10) lines++; }
    } catch (e) { error = errText(e); }
    add(`decompressionstream_cut_mate${mate}`, error === 'none' ? 'no_error' : 'error_at_end', `plain_bytes_before_end=${plain}; complete_records=${Math.floor(lines / 4)}; error=${error}`);
    let fplain = 0, flines = 0, ferror = 'none', partial = 'none';
    const { VerifyingGunzip } = await import(cfg.gzUrl);
    const g = new VerifyingGunzip(chunk => { fplain += chunk.length; for (const b of chunk) if (b === 10) flines++; });
    try { g.push(cut[mate].body, false); } catch (e) { partial = errText(e); }
    try { g.push(new Uint8Array(0), true); } catch (e) { ferror = errText(e); }
    add(`fflate_cut_mate${mate}`, ferror === 'none' ? 'no_error' : 'error_at_final', `plain_bytes=${fplain}; complete_records=${Math.floor(flines / 4)}; push_error=${partial}; final_error=${ferror}`);
    cut[mate].records = Math.floor(flines / 4);
  }
  // Index once for the kallisto tests.
  const txr = await fetch(cfg.txUrl);
  const tx = new Uint8Array(await txr.arrayBuffer());
  const ix = await s.runKallisto({ ...cfg.urls, args: s.KALLISTO_ARGS.index('t.fa', 't.idx'), files: [{ path: 't.fa', data: tx }], outputs: ['t.idx'], timeoutMs: cfg.timeoutMs });
  const idx = ix.files?.['t.idx'];
  add('index', idx ? 'pass' : 'fail', `exit=${ix.exitCode}; ms=${ix.ms.toFixed(0)}; error=${ix.error ?? 'none'}`);
  if (!idx) return rows;
  const quant = async (test, m1, m2, expect) => {
    const q = await s.runKallisto({ ...cfg.urls, args: s.KALLISTO_ARGS.quant('t.idx', 'kq', m1.name, m2.name),
      files: [{ path: 't.idx', data: idx }, { path: m1.name, data: m1.data }, { path: m2.name, data: m2.data }], outputs: ['kq/abundance.tsv', 'kq/run_info.json'], timeoutMs: cfg.timeoutMs });
    const info = q.files?.['kq/run_info.json'] ? s.runInfo(q.files['kq/run_info.json']) : {};
    const ab = q.files?.['kq/abundance.tsv'] ? await s.sha256(q.files['kq/abundance.tsv']) : 'none';
    const ok = !q.error && q.exitCode === 0 && (expect == null || info.n_processed === expect);
    add(test, ok ? 'pass' : 'fail', `exit=${q.exitCode}; error=${q.error ?? 'none'}; n_processed=${info.n_processed}; n_pseudoaligned=${info.n_pseudoaligned}; expected_pairs=${expect ?? 'n/a'}; abundance_sha256=${ab}; ms=${q.ms.toFixed(0)}; stderr_tail=${(q.stderr || '').slice(-600)}`);
    return { info, ab };
  };
  // kallisto on the raw cut gzip files, same byte count for both mates.
  await quant('kallisto_raw_cut_same_bytes', { name: 'c_1.fastq.gz', data: cut[1].body }, { name: 'c_2.fastq.gz', data: cut[2].body }, null);
  rows[rows.length - 1].detail += `; records_mate1=${cut[1].records}; records_mate2=${cut[2].records}`;
  // Different byte counts: mate 2 cut at 3/4 of the bytes.
  const short2 = cut[2].body.subarray(0, R2);
  let short2Records = 0;
  { const { VerifyingGunzip } = await import(cfg.gzUrl); const g = new VerifyingGunzip(c => { for (const b of c) if (b === 10) short2Records++; }); try { g.push(short2, false); } catch {} short2Records = Math.floor(short2Records / 4); }
  await quant('kallisto_raw_cut_different_bytes', { name: 'd_1.fastq.gz', data: cut[1].body }, { name: 'd_2.fastq.gz', data: short2 }, null);
  rows[rows.length - 1].detail += `; mate1_bytes=${cut[1].body.length}; mate2_bytes=${short2.length}; records_mate1=${cut[1].records}; records_mate2=${short2Records}`;
  // The N-pair method on Range streams: qmutfind keepFirstReads (fflate, packs gzip level 3) to N records per mate, then stop.
  const N = cfg.n;
  const viaRange = {};
  for (const mate of [1, 2]) {
    // The input is the bytes of the Range request above (the same bytes a live Range stream gives), read as a stream.
    const t = performance.now();
    const m = await s.prepMate(new Blob([cut[mate].body]).stream().getReader(), N, 'gz', `mate ${mate}`);
    viaRange[mate] = { ...m, ms: performance.now() - t, status: cut[mate].status };
  }
  const p1 = s.plainOf(viaRange[1].bytes), p2 = s.plainOf(viaRange[2].bytes);
  const pc = s.pairCheck(p1, p2);
  const plainSha = { 1: await s.sha256(p1), 2: await s.sha256(p2) };
  add('method_cut_bytes_keepFirstReads', viaRange[1].reads === N && viaRange[2].reads === N && pc.names_equal ? 'pass' : 'fail',
    `status=${viaRange[1].status}/${viaRange[2].status}; reads=${viaRange[1].reads}/${viaRange[2].reads}; cut=${viaRange[1].cut}/${viaRange[2].cut}; packed_bytes=${viaRange[1].bytes.length}/${viaRange[2].bytes.length}; ms=${viaRange[1].ms.toFixed(0)}/${viaRange[2].ms.toFixed(0)}; names_equal=${pc.names_equal}; first_diff=${pc.first_diff}; plain_sha256_1=${plainSha[1]}; plain_sha256_2=${plainSha[2]}`);
  await quant('kallisto_method_N_pairs', { name: 'n_1.fastq.gz', data: viaRange[1].bytes }, { name: 'n_2.fastq.gz', data: viaRange[2].bytes }, N);
  // Same method without a Range header: stream the whole URL, stop at N records, cancel.
  const full = {};
  try { for (const mate of [1, 2]) {
    const t = performance.now();
    let downloaded = 0;
    const r = await fetch(cfg.ena[mate], { mode: 'cors', credentials: 'omit' });
    const counting = r.body.pipeThrough(new TransformStream({ transform(chunk, ctl) { downloaded += chunk.length; ctl.enqueue(chunk); } }));
    const m = await s.prepMate(counting.getReader(), N, 'gz', `mate ${mate}`);
    full[mate] = { ...m, ms: performance.now() - t, status: r.status, downloaded, sha: await s.sha256(s.plainOf(m.bytes)) };
  } } catch (e) { add('method_stream_cancel_no_range', 'fail', errText(e)); full[1] = full[1] || {}; full[2] = full[2] || {}; }
  if (full[1].reads !== undefined && full[2].reads !== undefined) add('method_stream_cancel_no_range', full[1].reads === N && full[2].reads === N && full[1].sha === plainSha[1] && full[2].sha === plainSha[2] ? 'pass' : 'fail',
    `status=${full[1].status}/${full[2].status}; reads=${full[1].reads}/${full[2].reads}; bytes_read_before_cancel=${full[1].downloaded}/${full[2].downloaded}; ms=${full[1].ms.toFixed(0)}/${full[2].ms.toFixed(0)}; kB_per_s=${(full[1].downloaded / full[1].ms).toFixed(1)}/${(full[2].downloaded / full[2].ms).toFixed(1)}; plain_sha_equal_to_range_method=${full[1].sha === plainSha[1]}/${full[2].sha === plainSha[2]}`);
  // The same cut with DecompressionStream instead of fflate.
  const ds = {};
  for (const mate of [1, 2]) {
    try {
      const m = await s.keepFirstPlainDS(new Blob([cut[mate].body]).stream(), N);
      ds[mate] = { reads: m.reads, sha: await s.sha256(m.bytes) };
    } catch (e) { ds[mate] = { error: errText(e) }; }
  }
  add('method_cut_bytes_decompressionstream', ds[1].reads === N && ds[2].reads === N && ds[1].sha === plainSha[1] && ds[2].sha === plainSha[2] ? 'pass' : 'fail',
    `reads=${ds[1].reads}/${ds[2].reads}; plain_sha_equal_to_fflate=${ds[1].sha === plainSha[1]}/${ds[2].sha === plainSha[2]}; error=${ds[1].error ?? 'none'}/${ds[2].error ?? 'none'}`);
  // A range that is too short for the asked count (the stream ends first).
  const big = cut[1].records * 2;
  try {
    const m = await s.prepMate(new Blob([cut[1].body]).stream().getReader(), big, 'gz', 'mate 1');
    add('method_range_too_short', 'no_error', `asked=${big}; reads=${m.reads}; cut=${m.cut}`);
  } catch (e) { add('method_range_too_short', 'error', `asked=${big}; error=${errText(e)}`); }
  // The ENA filereport (read count) from the page.
  try {
    const r = await fetch(cfg.filereport, { mode: 'cors', credentials: 'omit' });
    add('ena_filereport', r.ok ? 'pass' : 'fail', `status=${r.status}; body=${(await r.text()).replace(/\n/g, ' | ').slice(0, 300)}`);
  } catch (e) { add('ena_filereport', 'fail', errText(e)); }
  return rows;
}

async function suiteTruncated(playwright) {
  const reps = Number(options.reps);
  const lines = ['engine\trep\ttest\tresult\tdetail\tbrowser_version'];
  const tsv = join(outDir, `truncated.${options.engine}.tsv`);
  for (let rep = 1; rep <= reps; rep++) {
    let session;
    try {
      session = await openBrowser(playwright);
      const port = server.address().port;
      const rows = await session.eval(truncatedInPage, { stepsUrl: `http://127.0.0.1:${port}/measure/steps.mjs`, gzUrl: `http://127.0.0.1:${port}/measure/qm/gz.js?v=__BUILD__`,
        urls: pageUrls(), timeoutMs: TIMEOUT_MS, txUrl: '/data/yeast.transcripts.fa', ena: ENA, n: Number(options['n-method']), rangeMiB: Number(options['range-mib']),
        filereport: `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,read_count,fastq_bytes&format=tsv` });
      for (const r of rows) lines.push([options.engine, rep, r.test, r.result, clean(r.detail), session.version].join('\t'));
    } catch (error) {
      lines.push([options.engine, rep, 'harness', errorOutcome(String(error)), clean(String(error?.message || error)), session?.version ?? ''].join('\t'));
    } finally { await closeBrowser(session); }
    await writeFile(tsv, lines.join('\n') + '\n');
    log(`truncated rep ${rep} done`);
  }
}

// ================================================================= suite genomes
const GENOME_COLS = ['engine', 'accession', 'rep', 'outcome', 'fna_gz_bytes', 'gff_gz_bytes', 'fna_bytes', 'gff_bytes',
  'gunzip_ds_s', 'gunzip_fflate_s', 'gunzip_equal', 'fna_members', 'gff_members', 'build_s', 'n_transcripts', 'n_genes', 'transcripts_bytes', 'transcripts_sha256',
  'index_s', 'index_heap_final_bytes', 'load_tree_pss_kb', 'build_tree_pss_peak_kb', 'index_tree_pss_peak_kb', 'node_rss_build_peak_kb', 'node_rss_index_peak_kb', 'node_maxrss_kb',
  'ncbi_fetch_s', 'ncbi_fetch_bytes', 'ncbi_sha_equal', 'browser_version', 'detail'];

async function genomeGunzipBuild(cfg) {
  const s = await import(cfg.stepsUrl);
  const { VerifyingGunzip } = await import(cfg.gzUrl);
  const out = {};
  const get = async url => { const r = await fetch(url); if (!r.ok) throw new Error(`fetch ${url}: HTTP ${r.status}`); return new Uint8Array(await r.arrayBuffer()); };
  const fnaGz = cfg.fnaGz || await get(cfg.fnaUrl), gffGz = cfg.gffGz || await get(cfg.gffUrl);
  out.fna_gz_bytes = fnaGz.length; out.gff_gz_bytes = gffGz.length;
  const fflate = gz => { const parts = []; const g = new VerifyingGunzip(c => parts.push(c)); g.push(gz, true); return { bytes: s.concat(parts), members: g.members }; };
  let t = performance.now();
  const fnaF = fflate(fnaGz), gffF = fflate(gffGz);
  out.gunzip_fflate_s = ((performance.now() - t) / 1000).toFixed(2);
  out.fna_members = fnaF.members; out.gff_members = gffF.members;
  let dsEqual = 'not_run';
  if (typeof DecompressionStream !== 'undefined') {
    t = performance.now();
    try {
      const ds = async gz => new Uint8Array(await new Response(new Blob([gz]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());
      const fnaD = await ds(fnaGz), gffD = await ds(gffGz);
      out.gunzip_ds_s = ((performance.now() - t) / 1000).toFixed(2);
      dsEqual = String((await s.sha256(fnaD)) === (await s.sha256(fnaF.bytes)) && (await s.sha256(gffD)) === (await s.sha256(gffF.bytes)));
    } catch (e) { dsEqual = `error: ${e?.name}: ${e?.message}`; }
  }
  out.gunzip_equal = dsEqual;
  out.fna_bytes = fnaF.bytes.length; out.gff_bytes = gffF.bytes.length;
  const wasm = cfg.isNode ? globalThis.measureWasm : await import(cfg.pkgUrl);
  if (!cfg.isNode) await wasm.default();
  t = performance.now();
  const gffText = new TextDecoder().decode(gffF.bytes);
  try {
    const tr = wasm.build_transcriptome(fnaF.bytes, gffText);
    out.build_s = ((performance.now() - t) / 1000).toFixed(2);
    const fasta = tr.fasta;
    out.n_transcripts = tr.n_transcripts; out.n_genes = tr.n_genes; out.transcripts_bytes = fasta.length;
    out.transcripts_sha256 = await s.sha256(fasta);
    tr.free?.();
    globalThis.measureTx = fasta;
  } catch (e) { out.build_s = ((performance.now() - t) / 1000).toFixed(2); out.error = `build_transcriptome: ${e?.message || e}`; }
  return out;
}
async function genomeIndex(cfg) {
  const s = await import(cfg.stepsUrl);
  const r = await s.runKallisto({ ...cfg.urls, args: s.KALLISTO_ARGS.index('g.fa', 'g.idx'), files: [{ path: 'g.fa', data: globalThis.measureTx }], outputs: ['g.idx'], timeoutMs: cfg.timeoutMs });
  globalThis.measureTx = null;
  return { ms: r.ms, exitCode: r.exitCode, error: r.error, missing: r.missing, heap: r.heap, stderrTail: (r.stderr || '').slice(-300) };
}
async function ncbiFetch(cfg) {
  const t = performance.now();
  const r = await fetch(cfg.url, { mode: 'cors', credentials: 'omit' });
  const b = new Uint8Array(await r.arrayBuffer());
  const digest = Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256', b)), x => x.toString(16).padStart(2, '0')).join('');
  return { status: r.status, bytes: b.length, ms: performance.now() - t, sha: digest };
}

function readGenomes() {
  return readFileSync(resolve(options.genomes), 'utf8').trim().split('\n').slice(1).map(l => {
    const [accession, fna, gff, fnaUrl, gffUrl, fnaSha, gffSha] = l.split('\t');
    return { accession, fna, gff, fnaUrl, gffUrl, fnaSha, gffSha };
  });
}

async function genomeChild() {
  const g = readGenomes().find(x => x.accession === options.genome);
  const row = { engine: 'node', accession: g.accession, rep: options.rep };
  let rssPeak = 0;
  const timer = setInterval(() => { rssPeak = Math.max(rssPeak, process.memoryUsage().rss); }, SAMPLE_MS);
  const phase = () => { const p = Math.max(rssPeak, process.memoryUsage().rss); rssPeak = 0; return kb(p); };
  try {
    const require = createRequire(import.meta.url);
    globalThis.measureWasm = require(join(pkgDir, 'qrnamap_wasm.js'));
    const cfg = { isNode: true, stepsUrl: pathToFileURL(join(here, 'steps.mjs')).href, gzUrl: pathToFileURL(join(here, 'qm', 'gz.js')).href + '?v=__BUILD__',
      urls: nodeUrls(), timeoutMs: TIMEOUT_MS,
      fnaGz: new Uint8Array(await readFile(join(dataDir, g.fna))), gffGz: new Uint8Array(await readFile(join(dataDir, g.gff))) };
    phase();
    Object.assign(row, await genomeGunzipBuild(cfg));
    row.node_rss_build_peak_kb = phase();
    if (row.error) { row.outcome = 'fail'; row.detail = row.error; }
    else {
      const ix = await genomeIndex(cfg);
      Object.assign(row, { index_s: s2(ix.ms), index_heap_final_bytes: ix.heap?.finalBytes, node_rss_index_peak_kb: phase() });
      row.outcome = ix.error || ix.exitCode !== 0 || ix.missing?.length ? errorOutcome(ix.error || 'fail') : 'finished';
      if (row.outcome !== 'finished') row.detail = ix.error || `exit ${ix.exitCode}; ${ix.stderrTail}`;
    }
  } catch (error) { row.detail = String(error?.stack || error).slice(0, 1200); row.outcome = errorOutcome(row.detail); }
  clearInterval(timer);
  row.node_maxrss_kb = process.resourceUsage().maxRSS;
  process.stdout.write('MEASURE_RESULT ' + JSON.stringify(row) + '\n');
}

async function genomeBrowserRun(playwright, g, rep) {
  const row = { engine: options.engine, accession: g.accession, rep };
  let session;
  try {
    session = await openBrowser(playwright);
    row.browser_version = session.version;
    const port = server.address().port;
    const cfg = { stepsUrl: `http://127.0.0.1:${port}/measure/steps.mjs`, gzUrl: `http://127.0.0.1:${port}/measure/qm/gz.js?v=__BUILD__`,
      pkgUrl: `http://127.0.0.1:${port}/pkg/qrnamap_wasm.js`, urls: pageUrls(), timeoutMs: TIMEOUT_MS, fnaUrl: `/data/${g.fna}`, gffUrl: `/data/${g.gff}` };
    if (rep === 1) {
      const a = await session.eval(ncbiFetch, { url: g.fnaUrl }), b = await session.eval(ncbiFetch, { url: g.gffUrl });
      Object.assign(row, { ncbi_fetch_s: s2(a.ms + b.ms), ncbi_fetch_bytes: a.bytes + b.bytes, ncbi_sha_equal: String(a.status === 200 && b.status === 200 && a.sha === g.fnaSha && b.sha === g.gffSha) });
    }
    await new Promise(done => setTimeout(done, 500));
    row.load_tree_pss_kb = sampleTree().pss;
    session.sampler.reset();
    Object.assign(row, await session.eval(genomeGunzipBuild, cfg));
    row.build_tree_pss_peak_kb = session.sampler.peak.pss;
    if (row.error) { row.outcome = 'fail'; row.detail = row.error; return row; }
    session.sampler.reset();
    const ix = await session.eval(genomeIndex, cfg);
    Object.assign(row, { index_s: s2(ix.ms), index_heap_final_bytes: ix.heap?.finalBytes, index_tree_pss_peak_kb: session.sampler.peak.pss });
    row.outcome = ix.error || ix.exitCode !== 0 || ix.missing?.length ? errorOutcome(ix.error || 'fail') : 'finished';
    if (row.outcome !== 'finished') row.detail = ix.error || `exit ${ix.exitCode}; ${ix.stderrTail}`;
  } catch (error) {
    row.detail = String(error?.message || error).slice(0, 1200); row.outcome = errorOutcome(row.detail);
  } finally { await closeBrowser(session); }
  return row;
}

async function suiteGenomes(playwright) {
  const reps = Number(options.reps);
  const tsv = join(outDir, `genomes.${options.engine}.tsv`);
  const lines = [GENOME_COLS.join('\t')];
  for (const g of readGenomes()) for (let rep = 1; rep <= reps; rep++) {
    log(`genome ${g.accession} rep ${rep} starting`);
    const row = options.engine === 'node'
      ? await runNodeChild(['--genomes', resolve(options.genomes), '--child', 'genome', '--genome', g.accession, '--rep', String(rep)], { engine: 'node', accession: g.accession, rep })
      : await genomeBrowserRun(playwright, g, rep);
    lines.push(toTsv(GENOME_COLS, row));
    await writeFile(tsv, lines.join('\n') + '\n');
    log(`genome ${g.accession} rep ${rep}: ${row.outcome} build_s=${row.build_s ?? ''} index_s=${row.index_s ?? ''} ${row.detail ?? ''}`);
  }
}

// ================================================================= suite determinism
// One browser (or Node process) per rep: one index build, then the first N pairs as gz and as plain, each quantified twice
// with that same index. The index and every abundance.tsv are written to OUT/determinism/ (artifact only, not git).
async function determinismInPage(cfg) {
  const s = await import(cfg.stepsUrl);
  const out = { quants: [] };
  const tx = globalThis.measureTx || new Uint8Array(await (await fetch(cfg.txUrl)).arrayBuffer());
  const ix = await s.runKallisto({ ...cfg.urls, args: s.KALLISTO_ARGS.index('yeast.transcripts.fa', 'yeast.idx'),
    files: [{ path: 'yeast.transcripts.fa', data: tx }], outputs: ['yeast.idx'], timeoutMs: cfg.timeoutMs });
  const idx = ix.files?.['yeast.idx'];
  if (!idx) { out.error = `index: ${ix.error || ix.exitCode}`; return out; }
  out.index_sha256 = await s.sha256(idx);
  out.index_b64 = cfg.keepIndex ? btoa(Array.from(idx, b => String.fromCharCode(b)).join('')) : '';
  for (const mode of ['gz', 'plain']) {
    const reads = {};
    for (const mate of [1, 2]) {
      const reader = cfg.open ? await cfg.open(mate) : (await fetch(cfg.readUrls[mate])).body.getReader();
      reads[mate] = (await s.prepMate(reader, cfg.pairs, mode, `mate ${mate}`)).bytes;
    }
    const ext = mode === 'gz' ? '.fastq.gz' : '.fastq';
    for (const k of [1, 2]) {
      const q = await s.runKallisto({ ...cfg.urls, args: s.KALLISTO_ARGS.quant('yeast.idx', 'kq', `r_1${ext}`, `r_2${ext}`),
        files: [{ path: 'yeast.idx', data: idx }, { path: `r_1${ext}`, data: reads[1] }, { path: `r_2${ext}`, data: reads[2] }],
        outputs: ['kq/abundance.tsv', 'kq/run_info.json'], timeoutMs: cfg.timeoutMs });
      const ab = q.files?.['kq/abundance.tsv'];
      out.quants.push({ mode, k, exitCode: q.exitCode, error: q.error, sha: ab ? await s.sha256(ab) : '', text: ab ? new TextDecoder().decode(ab) : '',
        input1_sha256: await s.sha256(reads[1]), ...(q.files?.['kq/run_info.json'] ? s.runInfo(q.files['kq/run_info.json']) : {}) });
    }
  }
  return out;
}

async function suiteDeterminism(playwright) {
  const reps = Number(options.reps), pairs = Number(options.counts.split(',')[0]);
  const dir = join(outDir, 'determinism');
  await mkdir(dir, { recursive: true });
  const lines = ['engine\tpairs\trep\tindex_sha256\tmode\tquant\texit_code\tn_processed\tn_pseudoaligned\tinput1_sha256\tabundance_sha256\tdetail'];
  for (let rep = 1; rep <= reps; rep++) {
    let session, result;
    try {
      if (options.engine === 'node') {
        globalThis.measureTx = new Uint8Array(await readFile(join(dataDir, 'yeast.transcripts.fa')));
        result = await determinismInPage({ stepsUrl: pathToFileURL(join(here, 'steps.mjs')).href, urls: nodeUrls(), timeoutMs: TIMEOUT_MS, pairs, keepIndex: true,
          open: async mate => Readable.toWeb(createReadStream(join(dataDir, `${RUN}_${mate}.fastq.gz`))).getReader() });
      } else {
        session = await openBrowser(playwright);
        const port = server.address().port;
        result = await session.eval(determinismInPage, { stepsUrl: `http://127.0.0.1:${port}/measure/steps.mjs`, urls: pageUrls(), timeoutMs: TIMEOUT_MS, pairs, keepIndex: true,
          txUrl: '/data/yeast.transcripts.fa', readUrls: { 1: `/data/${RUN}_1.fastq.gz`, 2: `/data/${RUN}_2.fastq.gz` } });
      }
    } catch (error) { result = { error: String(error?.message || error), quants: [] }; }
    finally { await closeBrowser(session); }
    if (result.index_b64) await writeFile(join(dir, `${options.engine}.r${rep}.yeast.idx`), Buffer.from(result.index_b64, 'base64'));
    if (result.error) lines.push([options.engine, pairs, rep, result.index_sha256 ?? '', '', '', '', '', '', '', '', clean(result.error)].join('\t'));
    for (const q of result.quants) {
      if (q.text) await writeFile(join(dir, `${options.engine}.r${rep}.${q.mode}.q${q.k}.abundance.tsv`), q.text);
      lines.push([options.engine, pairs, rep, result.index_sha256, q.mode, q.k, q.exitCode, q.n_processed, q.n_pseudoaligned, q.input1_sha256, q.sha, clean(q.error)].join('\t'));
    }
    await writeFile(join(outDir, `determinism.${options.engine}.tsv`), lines.join('\n') + '\n');
    log(`determinism rep ${rep} done: index ${result.index_sha256 ?? 'none'}`);
  }
}

// ================================================================= main
async function main() {
  if (options.child === 'quant') { await quantChild(); return 0; }
  if (options.child === 'genome') { await genomeChild(); return 0; }
  await mkdir(outDir, { recursive: true });
  const meta = { suite: options.suite, engine: options.engine, node: process.version, startedAt: new Date().toISOString(), sampleMs: SAMPLE_MS,
    timeoutMs: TIMEOUT_MS, guardKb: GUARD_KB, stopKb: STOP_KB, options, files: {} };
  for (const f of ['kallisto.mjs', 'kallisto.wasm', 'runtime/index.mjs', 'runtime/worker.mjs', 'runtime/execute.mjs']) meta.files[`site/vendor/kallisto/${f}`] = hash(await readFile(join(vendorDir, f)));
  for (const f of ['measure.mjs', 'steps.mjs', 'heap_wrap.mjs', 'qm/reads.js', 'qm/gz.js', 'qm/vendor/fflate.js']) meta.files[`measure/${f}`] = hash(await readFile(join(here, f)));
  if (pkgDir) for (const f of readdirSync(pkgDir)) meta.files[`pkg/${f}`] = hash(await readFile(join(pkgDir, f)));
  try { meta.memTotalKb = Number(/^MemTotal:\s+(\d+)/m.exec(readFileSync('/proc/meminfo', 'utf8'))[1]); meta.cpus = (await import('node:os')).cpus().length; } catch {}
  let playwright;
  if (options.engine !== 'node' || options.suite === 'truncated') {
    if (options.engine === 'node') throw new Error('suite truncated runs in browsers only');
    playwright = await import('playwright');
  }
  await startServer();
  try {
    if (options.suite === 'quant') await suiteQuant(playwright, meta);
    else if (options.suite === 'truncated') await suiteTruncated(playwright);
    else if (options.suite === 'genomes') await suiteGenomes(playwright);
    else if (options.suite === 'determinism') await suiteDeterminism(playwright);
    else throw new Error(`unknown suite ${options.suite}`);
  } finally {
    meta.finishedAt = new Date().toISOString();
    if (playwright) { try { const b = await playwright[options.engine].launch({ headless: true }); meta.browserVersion = b.version(); await b.close(); } catch {} }
    await writeFile(join(outDir, `${options.suite}.${options.engine}.meta.json`), JSON.stringify(meta, null, 2) + '\n');
    await new Promise(done => server.close(done));
  }
  return 0;
}

let status = 0;
try { status = await main(); } catch (error) { console.error(error?.stack || error); status = 2; }
process.exit(status);
