#!/usr/bin/env node
/**
 * Browser spike 0, Task 5 (spike 0c and 0d): cross-origin isolation and cross-origin downloads (CORS).
 * Runs inside the stajichlab/genomics-wasm fork (branch qrnamap-spike0), from the repository root.
 *
 * Usage:
 *   node spike0/isolation.mjs --engine chromium|firefox|webkit --data DIR --out DIR [--root DIR] [--modes none,coop-coep,service-worker]
 *
 * For each mode a new local server and a new browser launch are used:
 *   none            the local server sends no isolation header.
 *   coop-coep       Cross-Origin-Opener-Policy: same-origin, Cross-Origin-Embedder-Policy: require-corp,
 *                   Cross-Origin-Resource-Policy: same-origin (the headers of parity.mjs and of the upstream servers).
 *   service-worker  the local server sends no isolation header. The page registers a service worker (isolation-sw.js,
 *                   written after examples/rna-analysis/isolation-worker.js of the upstream repository). The worker adds the
 *                   three headers to every same-origin response. The page is then reloaded.
 * Per mode and browser:
 *   1. page facts: crossOriginIsolated, SharedArrayBuffer, shared WebAssembly.Memory, the same in a dedicated worker.
 *   2. tools: featureCounts -T 1, -T 2, -T 4 (pthreads), kallisto index and quant, HISAT2 build and align (same
 *      arguments as parity.mjs). A step that fails or hangs is a result: it is recorded with its error text.
 *   3. cross-origin fetches from the page (fetch, mode cors, no credentials) of the URLs in FETCHES. The page cannot
 *      read Access-Control-Allow-Origin or Cross-Origin-Resource-Policy from a response (they are not exposed to script),
 *      so the headers are taken from the Playwright network events ("seen"). The status and body come from the page.
 * Also, once per run: a probe with node:https (HEAD, or GET with Range) and an Origin header, no browser.
 * Failed, blocked and hung steps are results. Exit status: 0 when every step was attempted and recorded; 2 on a harness error.
 * Outputs (TSV): isolation_pages.tsv, isolation_tools.tsv, cors_results.tsv, cors_node_probe.tsv, meta-isolation.json.
 */
import { readFile, writeFile, mkdir } from 'node:fs/promises';
import { createServer } from 'node:http';
import { request as httpsRequest } from 'node:https';
import { createHash } from 'node:crypto';
import { resolve, join, relative, sep } from 'node:path';

const options = { engine: '', data: '', out: '', root: process.cwd(), modes: 'none,coop-coep,service-worker', 'tool-timeout-s': '240', 'fetch-timeout-s': '180' };
for (let i = 2; i < process.argv.length; i += 2) {
  const key = process.argv[i].replace(/^--/, '');
  if (!(key in options) || process.argv[i + 1] === undefined) throw new Error(`unknown or empty option ${process.argv[i]}`);
  options[key] = process.argv[i + 1];
}
if (!['chromium', 'firefox', 'webkit'].includes(options.engine)) throw new Error('--engine must be chromium, firefox or webkit');
if (!options.data || !options.out) throw new Error('--data and --out are required');
const root = resolve(options.root), dataDir = resolve(options.data), outDir = resolve(options.out);
const modes = options.modes.split(',').filter(Boolean);
const toolTimeoutMs = Number(options['tool-timeout-s']) * 1000;
const fetchTimeoutMs = Number(options['fetch-timeout-s']) * 1000;
const log = message => console.log(`[${new Date().toISOString()}] ${options.engine} ${message}`);
await mkdir(outDir, { recursive: true });

// URLs of the cross-origin tests (checked on 2026-11-21 with curl HEAD, see README).
const FETCHES = [
  { id: 'ncbi_fna', what: 'NCBI genome FASTA, S. cerevisiae R64, full GET', method: 'GET', expect: 3843460,
    url: 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz' },
  { id: 'ncbi_fna_range', what: 'NCBI genome FASTA, Range bytes=0-1023', method: 'GET', range: 'bytes=0-1023', expect: 1024,
    url: 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz' },
  { id: 'ncbi_gff', what: 'NCBI GFF3, S. cerevisiae R64, full GET', method: 'GET', expect: 2203839,
    url: 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz' },
  { id: 'ena_filereport', what: 'ENA portal API filereport of SRR3396381 (read_count), GET', method: 'GET',
    url: 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR3396381&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,read_count&format=tsv' },
  { id: 'ena_small_head', what: 'ENA FASTQ SRR22732385_1 (small), HEAD', method: 'HEAD',
    url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR227/085/SRR22732385/SRR22732385_1.fastq.gz' },
  { id: 'ena_small_get', what: 'ENA FASTQ SRR22732385_1 (small), full GET', method: 'GET', expect: 1280744,
    url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR227/085/SRR22732385/SRR22732385_1.fastq.gz' },
  { id: 'ena_large_head', what: 'ENA FASTQ SRR3396381_1 (669 MB), HEAD', method: 'HEAD',
    url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/001/SRR3396381/SRR3396381_1.fastq.gz' },
  { id: 'ena_large_range_start', what: 'ENA FASTQ SRR3396381_1, Range bytes=0-1048575', method: 'GET', range: 'bytes=0-1048575', expect: 1048576,
    url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/001/SRR3396381/SRR3396381_1.fastq.gz' },
  { id: 'ena_large_range_mid', what: 'ENA FASTQ SRR3396381_1, Range bytes=300000000-300065535', method: 'GET', range: 'bytes=300000000-300065535', expect: 65536,
    url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/001/SRR3396381/SRR3396381_1.fastq.gz' },
  { id: 'ensembl_gff3_mito', what: 'Ensembl GFF3, S. cerevisiae chromosome Mito (release 63 file name), full GET', method: 'GET',
    url: 'https://ftp.ensembl.org/pub/current_gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.63.chromosome.Mito.gff3.gz' },
  { id: 'ensembl_fa_mito', what: 'Ensembl FASTA, S. cerevisiae chromosome Mito, full GET', method: 'GET',
    url: 'https://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.Mito.fa.gz' },
  { id: 'ensembl_fa_checksums', what: 'Ensembl FASTA directory CHECKSUMS (3,741 bytes), full GET', method: 'GET', expect: 3741,
    url: 'https://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/dna/CHECKSUMS' },
  { id: 'ensemblgenomes_fungi', what: 'Ensembl Genomes fungi GFF3 directory listing (ftp.ensemblgenomes.org), GET', method: 'GET',
    url: 'https://ftp.ensemblgenomes.org/pub/fungi/current/gff3/' },
];

// Tool steps. inputs: [path, 'data:NAME' | 'kept:KEY']; keep: output files kept in the page for later steps.
const HT2 = Array.from({ length: 8 }, (_, i) => `yeast.hs.${i + 1}.ht2`);
const fcStep = threads => ({ tool: 'featurecounts', step: `count_T${threads}`, waitForExit: true,
  args: ['-T', String(threads), '-F', 'SAF', '-a', 'yeast.genes.saf', '-p', '--countReadPairs', '-s', '0', '-o', 'fc.tsv', 'yeast.sam'],
  inputs: [['yeast.genes.saf', 'data:yeast.genes.saf'], ['yeast.sam', 'data:native_yeast.sam']],
  outputs: ['fc.tsv', 'fc.tsv.summary'], keep: [], hash: { file: 'fc.tsv', mode: 'afterFirstLine' } });
const TOOL_STEPS = [
  fcStep(1), fcStep(2), fcStep(4),
  { tool: 'kallisto', step: 'index', waitForExit: false, args: ['index', '-i', 'yeast.idx', 'yeast.transcripts.fa'],
    inputs: [['yeast.transcripts.fa', 'data:yeast.transcripts.fa']], outputs: ['yeast.idx'], keep: ['yeast.idx'], hash: null },
  { tool: 'kallisto', step: 'quant', waitForExit: false, needs: 'kallisto index',
    args: ['quant', '-i', 'yeast.idx', '-o', 'kq', '-t', '1', '--plaintext', 'yeast_sub_1.fastq.gz', 'yeast_sub_2.fastq.gz'],
    inputs: [['yeast.idx', 'kept:yeast.idx'], ['yeast_sub_1.fastq.gz', 'data:yeast_sub_1.fastq.gz'], ['yeast_sub_2.fastq.gz', 'data:yeast_sub_2.fastq.gz']],
    outputs: ['kq/abundance.tsv', 'kq/run_info.json'], keep: [], hash: { file: 'kq/abundance.tsv', mode: 'all' } },
  { tool: 'hisat2', step: 'build', waitForExit: false, args: ['build', '-p', '1', 'yeast.fa', 'yeast.hs'],
    inputs: [['yeast.fa', 'data:yeast.fa']], outputs: HT2, keep: HT2, hash: null },
  { tool: 'hisat2', step: 'align', waitForExit: false, needs: 'hisat2 build',
    args: ['align', '-p', '1', '-x', 'yeast.hs', '-1', 'reads_1.fastq', '-2', 'reads_2.fastq', '-S', 'yeast.sam'],
    inputs: [...HT2.map(f => [f, `kept:${f}`]), ['reads_1.fastq', 'data:reads_1.fastq'], ['reads_2.fastq', 'data:reads_2.fastq']],
    outputs: ['yeast.sam'], keep: [], hash: { file: 'yeast.sam', mode: 'samBody' } },
];

const clean = value => String(value ?? '').replace(/[\t\r\n]+/g, ' ');
async function writeTsv(name, header, rows) {
  const text = [header.join('\t'), ...rows.map(row => header.map(key => clean(row[key])).join('\t'))].join('\n') + '\n';
  await writeFile(join(outDir, name), text);
}
const PAGE_HEADER = ['engine', 'browser_version', 'mode', 'run_date', 'phase', 'cross_origin_isolated', 'sab_defined', 'shared_memory_ok',
  'worker_cross_origin_isolated', 'worker_sab_defined', 'worker_shared_memory_ok', 'sw_controller', 'sw_note', 'error'];
const TOOL_HEADER = ['engine', 'browser_version', 'mode', 'cross_origin_isolated', 'tool', 'step', 'status', 'exit_code', 'seconds', 'args', 'output_hash', 'detail'];
const CORS_HEADER = ['engine', 'browser_version', 'mode', 'cross_origin_isolated', 'id', 'method', 'range', 'url', 'page_ok', 'page_status', 'response_type', 'page_content_length',
  'page_content_range', 'page_bytes_read', 'expected_bytes', 'bytes_match', 'body_magic', 'body_sha256', 'seen_status', 'seen_acao', 'seen_corp', 'seen_acah', 'seen_aceh', 'page_error', 'seen_failure'];
const PROBE_HEADER = ['id', 'method', 'range', 'url', 'status', 'acao', 'corp', 'acam', 'acah', 'aceh', 'content_length', 'content_range', 'content_type', 'accept_ranges', 'location', 'error'];
const pageRows = [], toolRows = [], corsRows = [], probeRows = [];
const flush = async () => {
  await writeTsv('isolation_pages.tsv', PAGE_HEADER, pageRows);
  await writeTsv('isolation_tools.tsv', TOOL_HEADER, toolRows);
  await writeTsv('cors_results.tsv', CORS_HEADER, corsRows);
  await writeTsv('cors_node_probe.tsv', PROBE_HEADER, probeRows);
};

// Node probe: headers as sent by each server for a request with an Origin header (no browser involved).
function probe(item) {
  return new Promise(done => {
    const row = { id: item.id, method: item.method, range: item.range || '', url: item.url };
    const headers = { Origin: 'http://127.0.0.1:8000', 'User-Agent': 'qrnamap-spike0-probe' };
    if (item.range) headers.Range = item.range;
    const req = httpsRequest(item.url, { method: item.method, headers, timeout: 60000 }, res => {
      const h = name => res.headers[name] ?? '';
      Object.assign(row, { status: res.statusCode, acao: h('access-control-allow-origin'), corp: h('cross-origin-resource-policy'),
        acam: h('access-control-allow-methods'), acah: h('access-control-allow-headers'), aceh: h('access-control-expose-headers'),
        content_length: h('content-length'), content_range: h('content-range'), content_type: h('content-type'),
        accept_ranges: h('accept-ranges'), location: h('location') });
      res.destroy(); // headers only
      done(row);
    });
    req.on('timeout', () => req.destroy(new Error('timeout after 60 s')));
    req.on('error', error => { row.error = `${error.code || ''} ${error.message}`.trim(); done(row); });
    req.end();
  });
}

const PROBE_WORKER = `self.onmessage = () => {
  let memory = false;
  try { new WebAssembly.Memory({ initial: 1, maximum: 2, shared: true }); memory = true; } catch (e) { memory = String(e.message || e); }
  postMessage({ coi: self.crossOriginIsolated, sab: typeof SharedArrayBuffer !== 'undefined', memory });
};`;
const ISOLATION_SW = `// Service worker that adds cross-origin isolation headers to same-origin responses.
// Written after examples/rna-analysis/isolation-worker.js of comparativechrono/genomics-wasm (same fetch logic and headers);
// the scope here is the whole local test origin, and there is no asset map.
self.addEventListener('install', event => event.waitUntil(self.skipWaiting()));
self.addEventListener('activate', event => event.waitUntil(self.clients.claim()));
self.addEventListener('fetch', event => {
  const request = event.request;
  const url = new URL(request.url);
  if (request.method !== 'GET' || url.origin !== self.location.origin) return;
  if (request.cache === 'only-if-cached' && request.mode !== 'same-origin') return;
  event.respondWith((async () => {
    const response = await fetch(request);
    if (response.status === 0) return response;
    const headers = new Headers(response.headers);
    headers.set('Cross-Origin-Opener-Policy', 'same-origin');
    headers.set('Cross-Origin-Embedder-Policy', 'require-corp');
    headers.set('Cross-Origin-Resource-Policy', 'same-origin');
    return new Response(response.body, { status: response.status, statusText: response.statusText, headers });
  })());
});`;

function inside(base, path) { return path === base || path.startsWith(base + sep); }
function makeServer(mode) {
  return createServer(async (request, response) => {
    if (mode === 'coop-coep') {
      response.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
      response.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
      response.setHeader('Cross-Origin-Resource-Policy', 'same-origin');
    }
    response.setHeader('Cache-Control', 'no-store');
    try {
      const name = decodeURIComponent(new URL(request.url, 'http://localhost').pathname);
      if (name === '/') { response.setHeader('Content-Type', 'text/html'); response.end('<!doctype html><meta charset="utf-8"><title>spike0 isolation</title>'); return; }
      if (name === '/probe-worker.js') { response.setHeader('Content-Type', 'text/javascript'); response.end(PROBE_WORKER); return; }
      if (name === '/isolation-sw.js') { response.setHeader('Content-Type', 'text/javascript'); response.setHeader('Service-Worker-Allowed', '/'); response.end(ISOLATION_SW); return; }
      if (name.startsWith('/data/')) {
        const path = resolve(dataDir, '.' + name.slice(5));
        if (!inside(dataDir, path)) { response.writeHead(403).end(); return; }
        response.setHeader('Content-Type', 'application/octet-stream');
        response.end(await readFile(path));
        return;
      }
      const path = resolve(root, '.' + name);
      if (!inside(root, path) || !['runtime/', 'dist/'].some(prefix => relative(root, path).startsWith(prefix))) { response.writeHead(404).end(); return; }
      response.setHeader('Content-Type', path.endsWith('.wasm') ? 'application/wasm' : 'text/javascript');
      response.end(await readFile(path));
    } catch (error) { log(`server error ${request.method} ${request.url}: ${error.message}`); response.writeHead(404).end(); }
  });
}

function withTimeout(promise, ms, text) {
  let timer;
  return Promise.race([promise, new Promise((_, reject) => { timer = setTimeout(() => reject(new Error(text)), ms); })]).finally(() => clearTimeout(timer));
}

// Facts of the current page and of a dedicated worker.
async function pageFacts(page) {
  return page.evaluate(async () => {
    const facts = { coi: globalThis.crossOriginIsolated, sab: typeof SharedArrayBuffer !== 'undefined', memory: false, worker: null };
    try { new WebAssembly.Memory({ initial: 1, maximum: 2, shared: true }); facts.memory = true; } catch (e) { facts.memory = String(e.message || e); }
    try {
      facts.worker = await new Promise((resolveWorker, rejectWorker) => {
        const worker = new Worker('/probe-worker.js');
        const timer = setTimeout(() => rejectWorker(new Error('worker timeout 20 s')), 20000);
        worker.onmessage = event => { clearTimeout(timer); worker.terminate(); resolveWorker(event.data); };
        worker.onerror = event => { clearTimeout(timer); rejectWorker(new Error(`worker error: ${event.message}`)); };
        worker.postMessage(1);
      });
    } catch (e) { facts.workerError = String(e.message || e); }
    return facts;
  });
}
function pushPage(base, phase, facts, extra = {}) {
  pageRows.push({ ...base, phase, cross_origin_isolated: facts.coi, sab_defined: facts.sab, shared_memory_ok: facts.memory,
    worker_cross_origin_isolated: facts.worker?.coi, worker_sab_defined: facts.worker?.sab, worker_shared_memory_ok: facts.worker?.memory,
    error: facts.workerError, ...extra });
}

// One tool step in the page. Inputs come from this server (/data/NAME) or from earlier steps kept in the page.
function pageStep(page, step) {
  const inputs = step.inputs.map(([path, src]) => ({ path, kind: src.slice(0, src.indexOf(':')), name: src.slice(src.indexOf(':') + 1) }));
  return page.evaluate(async ({ tool, args, inputs, outputs, keep, hash, waitForExit, timeoutMs }) => {
    const sha = async bytes => [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))].map(b => b.toString(16).padStart(2, '0')).join('');
    const { ToolRunner } = await import('/runtime/index.mjs');
    globalThis.__kept ??= new Map();
    const files = [];
    for (const input of inputs) {
      if (input.kind === 'kept') {
        const data = globalThis.__kept.get(input.name);
        if (!data) throw new Error(`kept file ${input.name} missing`);
        files.push({ path: input.path, data });
      } else {
        const response = await fetch(`/data/${input.name}`);
        if (!response.ok) throw new Error(`fetch /data/${input.name}: HTTP ${response.status}`);
        files.push({ path: input.path, data: new Uint8Array(await response.arrayBuffer()) });
      }
    }
    const runner = new ToolRunner({ moduleUrl: new URL(`/dist/${tool}/${tool}.mjs`, location.href), waitForExit });
    const result = await runner.run(args, { files, outputs, timeoutMs });
    const byName = new Map(result.files.map(file => [file.path.replace(/^\/work\//, ''), file.data]));
    for (const name of keep) if (byName.has(name)) globalThis.__kept.set(name, byName.get(name));
    let outputHash = '';
    if (hash && byName.has(hash.file)) {
      let bytes = byName.get(hash.file);
      if (hash.mode === 'afterFirstLine') bytes = bytes.subarray(bytes.indexOf(10) + 1); // line 1 holds the command line
      if (hash.mode === 'samBody') { let pos = 0; while (bytes[pos] === 64) pos = bytes.indexOf(10, pos) + 1; bytes = bytes.subarray(pos); } // drop @ header lines
      outputHash = `${hash.file}:${hash.mode}=${await sha(bytes)}`;
    }
    const tail = bytes => new TextDecoder().decode(bytes || new Uint8Array()).slice(-500);
    return { exitCode: result.exitCode, missing: result.missingOutputs, outputHash, stderrTail: tail(result.stderr) };
  }, { tool: step.tool, args: step.args, inputs, outputs: step.outputs, keep: step.keep, hash: step.hash, waitForExit: step.waitForExit, timeoutMs: toolTimeoutMs });
}

// One cross-origin fetch from the page.
function pageFetch(page, item) {
  return page.evaluate(async ({ url, method, range, timeoutMs }) => {
    const out = { ok: false };
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), timeoutMs);
    try {
      const response = await fetch(url, { method, headers: range ? { Range: range } : {}, mode: 'cors', credentials: 'omit', redirect: 'follow', signal: controller.signal });
      Object.assign(out, { status: response.status, type: response.type, contentLength: response.headers.get('content-length'), contentRange: response.headers.get('content-range') });
      if (method === 'GET') {
        const bytes = new Uint8Array(await response.arrayBuffer());
        out.bytes = bytes.length;
        out.magic = [...bytes.subarray(0, 2)].map(b => b.toString(16).padStart(2, '0')).join('');
        try { out.sha = [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))].map(b => b.toString(16).padStart(2, '0')).join(''); } catch {}
      } else out.bytes = 0;
      out.ok = true;
    } catch (error) { out.error = `${error.name}: ${error.message}`; }
    finally { clearTimeout(timer); }
    return out;
  }, { url: item.url, method: item.method, range: item.range || '', timeoutMs: fetchTimeoutMs });
}

async function runMode(playwright, mode) {
  log(`mode ${mode}: start`);
  const server = makeServer(mode);
  await new Promise(done => server.listen(0, '127.0.0.1', done));
  const origin = `http://127.0.0.1:${server.address().port}`;
  const browser = await playwright[options.engine].launch({ headless: true });
  const base = { engine: options.engine, browser_version: browser.version(), mode, run_date: new Date().toISOString() };
  try {
    const context = await browser.newContext();
    let page = await context.newPage();
    const seen = [], failed = [], pending = [];
    const attach = target => {
      target.on('console', message => log(`browser ${mode} ${message.type()}: ${message.text().slice(0, 300)}`));
      target.on('pageerror', error => log(`browser ${mode} pageerror: ${error.message}`));
      target.on('crash', () => log(`browser ${mode} page crashed`));
      target.on('response', response => pending.push((async () => {
        try { seen.push({ url: response.url(), method: response.request().method(), status: response.status(), headers: await response.allHeaders() }); }
        catch (error) { seen.push({ url: response.url(), method: response.request().method(), status: response.status(), headers: {}, note: error.message }); }
      })()));
      target.on('requestfailed', request => failed.push({ url: request.url(), method: request.method(), failure: request.failure()?.errorText || '' }));
    };
    attach(page);
    await page.goto(`${origin}/`);
    pushPage(base, 'first-load', await pageFacts(page));
    if (mode === 'service-worker') {
      let note = '';
      try {
        const reg = await withTimeout(page.evaluate(async () => {
          if (!('serviceWorker' in navigator)) return { error: 'navigator.serviceWorker is not available' };
          await navigator.serviceWorker.register('/isolation-sw.js', { scope: '/' });
          const ready = await navigator.serviceWorker.ready;
          return { scope: ready.scope, state: ready.active?.state, controller: !!navigator.serviceWorker.controller };
        }), 60000, 'service worker activation timeout 60 s');
        note = JSON.stringify(reg);
        pushPage(base, 'after-register-before-reload', await pageFacts(page), { sw_controller: reg.controller, sw_note: note });
        await page.reload();
        const controller = await page.evaluate(() => !!navigator.serviceWorker?.controller);
        pushPage(base, 'after-reload', await pageFacts(page), { sw_controller: controller, sw_note: note });
      } catch (error) {
        pushPage(base, 'service-worker-error', { coi: '', sab: '', memory: '' }, { sw_note: note, error: String(error.message || error) });
      }
    }
    const last = await pageFacts(page);
    const coi = last.coi;
    pushPage(base, 'tools-and-fetch', last);
    await flush();

    // Tools.
    const okSteps = new Set();
    for (const step of TOOL_STEPS) {
      const row = { ...base, cross_origin_isolated: coi, tool: step.tool, step: step.step, args: step.args.join(' ') };
      toolRows.push(row);
      if (step.needs && !okSteps.has(step.needs)) { Object.assign(row, { status: 'skipped', detail: `${step.needs} did not succeed` }); continue; }
      const started = Date.now();
      log(`mode ${mode}: ${step.tool} ${step.step} starting`);
      try {
        const result = await withTimeout(pageStep(page, step), toolTimeoutMs + 60000, `step exceeded ${toolTimeoutMs + 60000} ms (no result)`);
        row.seconds = ((Date.now() - started) / 1000).toFixed(2);
        row.exit_code = result.exitCode;
        row.output_hash = result.outputHash;
        if (result.exitCode !== 0 || result.missing.length) { row.status = 'fail'; row.detail = `exit ${result.exitCode}; missing ${result.missing.join(',') || 'none'}; stderr tail: ${result.stderrTail}`; }
        else { row.status = 'ok'; okSteps.add(`${step.tool} ${step.step}`); }
      } catch (error) {
        row.seconds = ((Date.now() - started) / 1000).toFixed(2);
        Object.assign(row, { status: 'error', detail: String(error?.message || error).slice(0, 800) });
        if (page.isClosed()) { log('page closed; a new page is opened'); page = await context.newPage(); attach(page); await page.goto(`${origin}/`); }
      }
      log(`mode ${mode}: ${step.tool} ${step.step} ${row.status} after ${row.seconds} s ${row.status === 'ok' ? '' : row.detail}`);
      await flush();
    }

    // Cross-origin fetches.
    for (const item of FETCHES) {
      seen.length = 0; failed.length = 0; pending.length = 0;
      const row = { ...base, cross_origin_isolated: coi, id: item.id, method: item.method, range: item.range || '', url: item.url, expected_bytes: item.expect ?? '' };
      corsRows.push(row);
      log(`mode ${mode}: fetch ${item.id}`);
      let result;
      try { result = await withTimeout(pageFetch(page, item), fetchTimeoutMs + 30000, 'fetch exceeded its time limit (no result)'); }
      catch (error) { result = { ok: false, error: String(error.message || error) }; }
      await Promise.all(pending);
      const sameUrl = event => event.url.split('#')[0] === item.url && event.method === item.method;
      const event = seen.find(sameUrl) || seen.find(e => e.url.split('#')[0] === item.url);
      const failure = failed.find(f => f.url.split('#')[0] === item.url);
      Object.assign(row, { page_ok: result.ok, page_status: result.status, response_type: result.type, page_content_length: result.contentLength,
        page_content_range: result.contentRange, page_bytes_read: result.bytes, body_magic: result.magic, body_sha256: result.sha, page_error: result.error,
        bytes_match: item.expect === undefined || result.bytes === undefined ? '' : result.bytes === item.expect,
        seen_status: event?.status, seen_acao: event?.headers['access-control-allow-origin'], seen_corp: event?.headers['cross-origin-resource-policy'],
        seen_acah: event?.headers['access-control-allow-headers'], seen_aceh: event?.headers['access-control-expose-headers'],
        seen_failure: failure?.failure });
      await flush();
    }
  } finally {
    await browser.close().catch(() => {});
    await new Promise(done => server.close(done));
    log(`mode ${mode}: done`);
  }
}

let harnessError = false;
const meta = { engine: options.engine, node: process.version, startedAt: new Date().toISOString(), modes, fetches: FETCHES.map(f => ({ id: f.id, url: f.url, method: f.method, range: f.range || '' })) };
try {
  meta.runtimeSha256 = {};
  for (const name of ['index.mjs', 'worker.mjs', 'execute.mjs']) meta.runtimeSha256[name] = createHash('sha256').update(await readFile(join(root, 'runtime', name))).digest('hex');
  for (const item of FETCHES) { const row = await probe(item); probeRows.push(row); log(`node probe ${item.id}: status ${row.status ?? ''} acao ${row.acao ?? ''} ${row.error ?? ''}`); }
  await flush();
  const playwright = await import('playwright');
  for (const mode of modes) {
    try { await runMode(playwright, mode); }
    catch (error) { log(`mode ${mode}: harness error ${error?.stack || error}`); pageRows.push({ engine: options.engine, mode, phase: 'mode-error', error: String(error?.message || error).slice(0, 500) }); await flush(); }
  }
} catch (error) {
  harnessError = true;
  meta.harnessError = String(error?.stack || error);
  console.error(error?.stack || error);
} finally {
  meta.finishedAt = new Date().toISOString();
  await flush();
  await writeFile(join(outDir, 'meta-isolation.json'), JSON.stringify(meta, null, 2) + '\n');
}
process.exit(harnessError ? 2 : 0);
