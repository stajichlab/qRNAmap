#!/usr/bin/env node
/**
 * Browser spike 0, Task 4: index build time and peak memory of the WASM kallisto and HISAT2 packages,
 * in Node or in a Playwright browser. Runs inside the stajichlab/genomics-wasm fork (branch qrnamap-spike0),
 * from the repository root, with the upstream runtime (runtime/index.mjs ToolRunner) and dist/<tool>/.
 *
 * Usage:
 *   node spike0/index_size.mjs --engine node|chromium|firefox|webkit --data DIR --out DIR
 *        [--cases hisat2_yeast,hisat2_afum,kallisto_afum,kallisto_yeast] [--reps 3] [--root DIR]
 *
 * Each case runs --reps times. Every run uses a fresh process: a child Node process (engine node) or a
 * fresh browser launch (browsers). So the peak values of one run do not carry into the next run.
 *
 * Memory probes (all values are written; none is copied from another engine):
 *   wasm_heap_*      Module.HEAPU8.length before and after callMain, read in the tool worker by
 *                    index_size_wrap.mjs. Linear memory never shrinks, so "final" is the peak linear memory.
 *                    It does not include MEMFS files (JS typed arrays) or engine overhead.
 *   node_maxrss_kb   engine node: process.resourceUsage().maxRSS of the child process (kernel peak RSS of the
 *                    whole process: main thread, tool worker, WASM memory, input and output copies).
 *   node_rss_peak_kb engine node: peak of process.memoryUsage().rss, sampled every SAMPLE_MS.
 *   tree_pss_*_kb    browsers: sum of Pss (/proc/PID/smaps_rollup) over all processes that descend from this
 *                    Node process (the browser and its child processes), sampled every SAMPLE_MS.
 *                    baseline = after the page load, before the tool run; peak = highest sample during the run.
 *                    A peak shorter than SAMPLE_MS can be missed.
 *   tree_rss_*_kb    same processes: sum of VmRSS (shared pages counted once per process) and the largest
 *                    single-process VmRSS.
 *   cdp_js_heap_*    Chromium only: CDP Performance.getMetrics JSHeapUsedSize and JSHeapTotalSize, sampled every
 *                    second. This covers the page main thread only, not the tool worker and not WASM memory.
 * Times: run_s = duration of ToolRunner.run (worker start, module compile, input copy, callMain, output copy);
 *        callmain_s = duration of callMain alone (from the wrapper).
 * Outcome: finished | fail (exit code not 0 or output missing) | oom (error text names memory) | error |
 *          crash (page or process crash) | timeout. A failed run is a result; the next run still starts.
 * Exit status: 0 when every planned run was attempted and recorded; 2 on a harness error.
 */
import { readFile, writeFile, mkdir } from 'node:fs/promises';
import { readFileSync, readdirSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { createServer } from 'node:http';
import { spawn } from 'node:child_process';
import { resolve, join, relative, sep, dirname } from 'node:path';
import { pathToFileURL, fileURLToPath } from 'node:url';

const options = { engine: 'node', cases: 'hisat2_yeast,hisat2_afum,kallisto_afum,kallisto_yeast', reps: '3', data: '', out: '', root: process.cwd(), child: '', rep: '' };
for (let i = 2; i < process.argv.length; i += 2) {
  const key = process.argv[i].replace(/^--/, '');
  if (!(key in options) || process.argv[i + 1] === undefined) throw new Error(`unknown or empty option ${process.argv[i]}`);
  options[key] = process.argv[i + 1];
}
const ENGINES = ['node', 'chromium', 'firefox', 'webkit'];
if (!ENGINES.includes(options.engine)) throw new Error(`engine must be one of ${ENGINES}`);
if (!options.data || !options.out) throw new Error('--data and --out are required');
const root = resolve(options.root), dataDir = resolve(options.data), outDir = resolve(options.out);
const SAMPLE_MS = 250;
const TIMEOUT_MS = 45 * 60 * 1000; // per run; native builds take 9 to 25 s (native_times.tsv)
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const log = message => console.error(`[${new Date().toISOString()}] ${options.engine} ${message}`);
const wrapFile = join(dirname(fileURLToPath(import.meta.url)), 'index_size_wrap.mjs');

// Arguments equal the native arguments (native_refs.sh), apart from file paths.
const HT2 = prefix => Array.from({ length: 8 }, (_, i) => `${prefix}.${i + 1}.ht2`);
const CASES = {
  hisat2_yeast: { tool: 'hisat2', input: 'yeast.fa', args: ['build', '-p', '1', 'yeast.fa', 'yeast.hs'], outputs: HT2('yeast.hs'), native: 'hisat2_build_yeast' },
  hisat2_afum: { tool: 'hisat2', input: 'afum.fa', args: ['build', '-p', '1', 'afum.fa', 'afum.hs'], outputs: HT2('afum.hs'), native: 'hisat2_build_afum' },
  kallisto_afum: { tool: 'kallisto', input: 'afum.transcripts.fa', args: ['index', '-i', 'afum.idx', 'afum.transcripts.fa'], outputs: ['afum.idx'], native: '' },
  kallisto_yeast: { tool: 'kallisto', input: 'yeast.transcripts.fa', args: ['index', '-i', 'yeast.idx', 'yeast.transcripts.fa'], outputs: ['yeast.idx'], native: 'kallisto_index' },
};
const COLUMNS = ['engine', 'case', 'tool', 'input', 'input_bytes', 'rep', 'outcome', 'exit_code', 'run_s', 'callmain_s',
  'wasm_heap_initial_bytes', 'wasm_heap_final_bytes', 'node_maxrss_kb', 'node_rss_peak_kb',
  'tree_pss_baseline_kb', 'tree_pss_peak_kb', 'tree_rss_sum_peak_kb', 'tree_rss_max_proc_kb', 'tree_procs_peak',
  'cdp_js_heap_used_peak_bytes', 'cdp_js_heap_total_peak_bytes', 'samples', 'index_bytes', 'outputs', 'browser_version', 'detail'];

function classify(row) {
  if (row.outcome) return row.outcome;
  if (row.exit_code === 0 && !row.missing) return 'finished';
  return 'fail';
}
function errorOutcome(text) {
  if (/\bOOM\b|out of memory|Cannot enlarge memory|bad_alloc|Array buffer allocation failed|RangeError: .*memory/i.test(text)) return 'oom';
  if (/crash|Target (page, context or browser )?closed|has been closed|disconnected/i.test(text)) return 'crash';
  if (/exceeded \d+ ms|TimeoutError/.test(text)) return 'timeout';
  return 'error';
}
const kb = bytes => Math.round(bytes / 1024);

// ---------- child mode (engine node): one run in this process ----------
async function childRun() {
  const spec = CASES[options.child];
  const row = { engine: 'node', case: options.child, tool: spec.tool, input: spec.input, rep: options.rep };
  let rssPeak = process.memoryUsage().rss, samples = 0;
  const timer = setInterval(() => { samples++; rssPeak = Math.max(rssPeak, process.memoryUsage().rss); }, SAMPLE_MS);
  let broadcast = null;
  const channel = new BroadcastChannel('spike0-heap');
  channel.onmessage = event => { broadcast = event.data; };
  try {
    const { ToolRunner } = await import(pathToFileURL(join(root, 'runtime', 'index.mjs')).href);
    const data = await readFile(join(dataDir, spec.input));
    row.input_bytes = data.length;
    const toolUrl = pathToFileURL(join(root, 'dist', spec.tool, `${spec.tool}.mjs`)).href;
    const moduleUrl = new URL(`${pathToFileURL(wrapFile).href}?tool=${encodeURIComponent(toolUrl)}`);
    const runner = new ToolRunner({ moduleUrl, waitForExit: false });
    const started = performance.now();
    try {
      const result = await runner.run(spec.args, { files: [{ path: spec.input, data }], outputs: [...spec.outputs, 'spike0_heap.json'], timeoutMs: TIMEOUT_MS });
      row.run_s = ((performance.now() - started) / 1000).toFixed(2);
      row.exit_code = result.exitCode;
      const missing = result.missingOutputs.filter(p => !p.endsWith('spike0_heap.json'));
      row.missing = missing.length > 0;
      let heap = null, total = 0;
      const outputs = [];
      for (const file of result.files) {
        const name = file.path.replace(/^\/work\//, '');
        if (name === 'spike0_heap.json') { heap = JSON.parse(new TextDecoder().decode(file.data)); continue; }
        total += file.data.length;
        outputs.push(`${name}=${file.data.length}:${hash(file.data)}`);
      }
      row.index_bytes = total;
      row.outputs = outputs.join(';');
      row.heap = heap;
      await mkdir(join(outDir, 'logs'), { recursive: true });
      await writeFile(join(outDir, 'logs', `${options.child}.r${options.rep}.stderr.txt`), result.stderr);
      if (row.exit_code !== 0 || row.missing) row.detail = `exit ${row.exit_code}; missing ${missing.join(',') || 'none'}; stderr tail: ${new TextDecoder().decode(result.stderr).slice(-400)}`;
    } catch (error) {
      row.run_s = ((performance.now() - started) / 1000).toFixed(2);
      row.detail = String(error?.stack || error).slice(0, 1200);
      row.outcome = errorOutcome(row.detail);
    }
  } finally {
    await new Promise(done => setTimeout(done, 200)); // let a last BroadcastChannel message arrive
    clearInterval(timer);
    channel.close();
  }
  row.heap = row.heap || broadcast;
  row.node_rss_peak_kb = kb(rssPeak);
  row.node_maxrss_kb = process.resourceUsage().maxRSS; // kilobytes
  row.samples = samples;
  process.stdout.write('SPIKE0_RESULT ' + JSON.stringify(row) + '\n');
}

function runNodeChild(caseName, rep) {
  return new Promise(done => {
    const args = [fileURLToPath(import.meta.url), '--engine', 'node', '--child', caseName, '--rep', String(rep), '--data', dataDir, '--out', outDir, '--root', root];
    const child = spawn(process.execPath, args, { stdio: ['ignore', 'pipe', 'pipe'] });
    let stdout = '', stderr = '';
    child.stdout.on('data', chunk => { stdout += chunk; });
    child.stderr.on('data', chunk => { stderr += chunk; process.stderr.write(chunk); });
    child.on('close', (code, signal) => {
      const line = stdout.split('\n').find(l => l.startsWith('SPIKE0_RESULT '));
      if (line) { done(JSON.parse(line.slice(14))); return; }
      const detail = `child exit ${code} signal ${signal}; stderr tail: ${stderr.slice(-800)}`;
      done({ engine: 'node', case: caseName, tool: CASES[caseName].tool, input: CASES[caseName].input, rep, outcome: signal || code === 134 || code === 137 ? 'crash' : errorOutcome(detail), detail });
    });
  });
}

// ---------- browser mode ----------
function procTree() {
  // Processes that descend from this Node process (not this process itself).
  const parent = new Map();
  for (const name of readdirSyncSafe('/proc')) {
    if (!/^\d+$/.test(name)) continue;
    try {
      const stat = readFileSync(`/proc/${name}/stat`, 'utf8');
      const fields = stat.slice(stat.lastIndexOf(')') + 2).split(' ');
      parent.set(Number(name), Number(fields[1]));
    } catch {}
  }
  const result = [];
  for (const pid of parent.keys()) {
    let p = parent.get(pid), depth = 0;
    while (p && p !== process.pid && depth++ < 64) p = parent.get(p);
    if (p === process.pid) result.push(pid);
  }
  return result;
}
function readdirSyncSafe(path) { try { return readdirSync(path); } catch { return []; } }
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

let server;
function startServer() {
  server = createServer(async (request, response) => {
    response.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
    response.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
    response.setHeader('Cross-Origin-Resource-Policy', 'same-origin');
    try {
      const name = decodeURIComponent(new URL(request.url, 'http://localhost').pathname);
      if (name === '/') { response.setHeader('Content-Type', 'text/html'); response.end('<!doctype html><meta charset="utf-8"><title>spike0 index size</title>'); return; }
      if (name === '/spike0/index_size_wrap.mjs') { response.setHeader('Content-Type', 'text/javascript'); response.end(await readFile(wrapFile)); return; }
      if (name.startsWith('/data/')) {
        const path = resolve(dataDir, '.' + name.slice(5));
        if (path !== dataDir && !path.startsWith(dataDir + sep)) { response.writeHead(403).end(); return; }
        response.setHeader('Content-Type', 'application/octet-stream');
        response.end(await readFile(path));
        return;
      }
      const path = resolve(root, '.' + name);
      if (!path.startsWith(root + sep) || !['runtime/', 'dist/'].some(prefix => relative(root, path).startsWith(prefix))) { response.writeHead(404).end(); return; }
      response.setHeader('Content-Type', path.endsWith('.wasm') ? 'application/wasm' : 'text/javascript');
      response.end(await readFile(path));
    } catch (error) { log(`server error ${request.method} ${request.url}: ${error.message}`); response.writeHead(404).end(); }
  });
  return new Promise(done => server.listen(0, '127.0.0.1', done));
}

async function browserRun(playwright, caseName, rep) {
  const spec = CASES[caseName];
  const row = { engine: options.engine, case: caseName, tool: spec.tool, input: spec.input, rep };
  row.input_bytes = (await readFile(join(dataDir, spec.input))).length;
  let browser, cdp, treeTimer, cdpTimer;
  const peak = { pss: 0, rss: 0, maxRss: 0, procs: 0, pssMissing: 0, samples: 0, jsUsed: 0, jsTotal: 0 };
  try {
    browser = await playwright[options.engine].launch({ headless: true });
    row.browser_version = browser.version();
    const page = await browser.newPage();
    page.on('console', message => log(`browser ${message.type()}: ${message.text()}`));
    page.on('pageerror', error => log(`browser pageerror: ${error.message}`));
    let crashed;
    const crash = new Promise((_, reject) => { crashed = reject; });
    crash.catch(() => {});
    page.on('crash', () => { log('browser page crashed'); crashed(new Error('page crash event')); });
    browser.on('disconnected', () => crashed(new Error('browser disconnected')));
    await page.goto(`http://127.0.0.1:${server.address().port}/`);
    const isolated = await page.evaluate(() => crossOriginIsolated);
    // Load the inputs into the page first, so the baseline includes the input copy in the page.
    await page.evaluate(async input => {
      const response = await fetch(`/data/${input}`);
      if (!response.ok) throw new Error(`fetch ${input}: HTTP ${response.status}`);
      globalThis.spike0Input = new Uint8Array(await response.arrayBuffer());
    }, spec.input);
    await new Promise(done => setTimeout(done, 1000));
    const base = sampleTree();
    row.tree_pss_baseline_kb = base.pss;
    if (options.engine === 'chromium') {
      cdp = await page.context().newCDPSession(page);
      await cdp.send('Performance.enable');
      const readCdp = async () => {
        try {
          const { metrics } = await cdp.send('Performance.getMetrics');
          const get = name => metrics.find(m => m.name === name)?.value ?? 0;
          peak.jsUsed = Math.max(peak.jsUsed, get('JSHeapUsedSize'));
          peak.jsTotal = Math.max(peak.jsTotal, get('JSHeapTotalSize'));
        } catch {}
      };
      await readCdp();
      cdpTimer = setInterval(readCdp, 1000);
    }
    const sample = () => {
      const s = sampleTree();
      peak.samples++;
      peak.pss = Math.max(peak.pss, s.pss); peak.rss = Math.max(peak.rss, s.rss);
      peak.maxRss = Math.max(peak.maxRss, s.maxRss); peak.procs = Math.max(peak.procs, s.procs);
      peak.pssMissing = Math.max(peak.pssMissing, s.pssMissing);
    };
    sample();
    treeTimer = setInterval(sample, SAMPLE_MS);
    const evaluation = page.evaluate(async ({ tool, input, args, outputs, timeoutMs }) => {
      const { ToolRunner } = await import('/runtime/index.mjs');
      let broadcast = null;
      const channel = new BroadcastChannel('spike0-heap');
      channel.onmessage = event => { broadcast = event.data; };
      const toolUrl = new URL(`/dist/${tool}/${tool}.mjs`, location.href).href;
      const moduleUrl = new URL(`/spike0/index_size_wrap.mjs?tool=${encodeURIComponent(toolUrl)}`, location.href);
      const runner = new ToolRunner({ moduleUrl, waitForExit: false });
      const started = performance.now();
      const result = { };
      try {
        const run = await runner.run(args, { files: [{ path: input, data: globalThis.spike0Input }], outputs: [...outputs, 'spike0_heap.json'], timeoutMs });
        result.runMs = performance.now() - started;
        result.exitCode = run.exitCode;
        result.missing = run.missingOutputs.filter(p => !p.endsWith('spike0_heap.json'));
        result.outputs = [];
        result.indexBytes = 0;
        for (const file of run.files) {
          const name = file.path.replace(/^\/work\//, '');
          if (name === 'spike0_heap.json') { result.heap = JSON.parse(new TextDecoder().decode(file.data)); continue; }
          const digest = Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256', file.data)), b => b.toString(16).padStart(2, '0')).join('');
          result.indexBytes += file.data.length;
          result.outputs.push(`${name}=${file.data.length}:${digest}`);
        }
        result.stderr = new TextDecoder().decode(run.stderr);
      } catch (error) {
        result.runMs = performance.now() - started;
        result.error = String(error?.stack || error).slice(0, 1200);
      }
      await new Promise(done => setTimeout(done, 200));
      channel.close();
      result.heap = result.heap || broadcast;
      return result;
    }, { tool: spec.tool, input: spec.input, args: spec.args, outputs: spec.outputs, timeoutMs: TIMEOUT_MS });
    let guard;
    const result = await Promise.race([evaluation, crash, new Promise((_, reject) => {
      guard = setTimeout(() => reject(new Error(`browser run exceeded ${TIMEOUT_MS + 120000} ms`)), TIMEOUT_MS + 120000);
    })]).finally(() => clearTimeout(guard));
    row.run_s = (result.runMs / 1000).toFixed(2);
    row.heap = result.heap;
    if (result.error) { row.detail = result.error; row.outcome = errorOutcome(result.error); }
    else {
      row.exit_code = result.exitCode;
      row.missing = result.missing.length > 0;
      row.index_bytes = result.indexBytes;
      row.outputs = result.outputs.join(';');
      await mkdir(join(outDir, 'logs'), { recursive: true });
      await writeFile(join(outDir, 'logs', `${caseName}.r${rep}.stderr.txt`), result.stderr);
      if (row.exit_code !== 0 || row.missing) row.detail = `exit ${row.exit_code}; missing ${result.missing.join(',') || 'none'}; stderr tail: ${result.stderr.slice(-400)}`;
    }
    row.isolated = isolated;
  } catch (error) {
    row.detail = String(error?.stack || error).slice(0, 1200);
    row.outcome = errorOutcome(row.detail);
  } finally {
    clearInterval(treeTimer); clearInterval(cdpTimer);
    row.tree_pss_peak_kb = peak.pss; row.tree_rss_sum_peak_kb = peak.rss; row.tree_rss_max_proc_kb = peak.maxRss;
    row.tree_procs_peak = peak.procs; row.samples = peak.samples;
    if (peak.pssMissing) row.detail = `${row.detail ? row.detail + '; ' : ''}Pss not readable for ${peak.pssMissing} processes (Rss used)`;
    if (options.engine === 'chromium') { row.cdp_js_heap_used_peak_bytes = peak.jsUsed; row.cdp_js_heap_total_peak_bytes = peak.jsTotal; }
    if (browser) await browser.close().catch(() => {});
  }
  return row;
}

function toTsvRow(row) {
  const flat = { ...row, outcome: classify(row) };
  if (row.heap) {
    flat.wasm_heap_initial_bytes = row.heap.initialBytes; flat.wasm_heap_final_bytes = row.heap.finalBytes;
    flat.callmain_s = row.heap.callMainMs != null ? (row.heap.callMainMs / 1000).toFixed(2) : '';
    if (row.heap.threw) flat.detail = `${flat.detail ? flat.detail + '; ' : ''}callMain threw: ${row.heap.threw}`;
  }
  const clean = value => String(value ?? '').replace(/[\t\r\n]+/g, ' ');
  return COLUMNS.map(key => clean(flat[key])).join('\t');
}

async function main() {
  if (options.child) { await childRun(); return 0; }
  await mkdir(outDir, { recursive: true });
  const cases = options.cases.split(',').filter(Boolean);
  for (const c of cases) if (!CASES[c]) throw new Error(`unknown case ${c}`);
  const reps = Number(options.reps);
  const meta = { engine: options.engine, node: process.version, startedAt: new Date().toISOString(), sampleMs: SAMPLE_MS, timeoutMs: TIMEOUT_MS, cases, reps, builds: {}, inputs: {}, wrapper: hash(await readFile(wrapFile)) };
  for (const tool of new Set(cases.map(c => CASES[c].tool))) {
    meta.builds[tool] = {};
    for (const name of [`${tool}.mjs`, `${tool}.wasm`]) meta.builds[tool][name] = hash(await readFile(join(root, 'dist', tool, name)));
  }
  for (const c of cases) { const bytes = await readFile(join(dataDir, CASES[c].input)); meta.inputs[CASES[c].input] = { bytes: bytes.length, sha256: hash(bytes) }; }
  meta.runtime = {};
  for (const name of ['index.mjs', 'worker.mjs', 'execute.mjs']) meta.runtime[name] = hash(await readFile(join(root, 'runtime', name)));
  try { meta.memTotalKb = Number(/^MemTotal:\s+(\d+)/m.exec(readFileSync('/proc/meminfo', 'utf8'))[1]); meta.cpus = (await import('node:os')).cpus().length; } catch {}
  const tsvPath = join(outDir, `index_size.${options.engine}.tsv`);
  const lines = [COLUMNS.join('\t')];
  const raw = [];
  let playwright;
  if (options.engine !== 'node') { playwright = await import('playwright'); await startServer(); }
  try {
    for (const c of cases) for (let rep = 1; rep <= reps; rep++) {
      log(`${c} rep ${rep} starting: ${CASES[c].tool} ${CASES[c].args.join(' ')}`);
      const row = options.engine === 'node' ? await runNodeChild(c, rep) : await browserRun(playwright, c, rep);
      row.engine = options.engine;
      raw.push(row);
      lines.push(toTsvRow(row));
      await writeFile(tsvPath, lines.join('\n') + '\n');
      log(`${c} rep ${rep} ${classify(row)} run_s=${row.run_s ?? ''} heap=${row.heap?.finalBytes ?? ''} ${row.detail ?? ''}`);
    }
  } finally {
    meta.finishedAt = new Date().toISOString();
    await writeFile(join(outDir, `index_size.${options.engine}.meta.json`), JSON.stringify(meta, null, 2) + '\n');
    await writeFile(join(outDir, `index_size.${options.engine}.raw.json`), JSON.stringify(raw, null, 2) + '\n');
    if (server) await new Promise(done => server.close(done));
  }
  return 0;
}

let status = 0;
try { status = await main(); } catch (error) { console.error(error?.stack || error); status = 2; }
process.exit(status);
