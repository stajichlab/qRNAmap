// Peak memory and time of the real genome step: the page downloads Fulvia fulva Race5_Kim (GCF_020509005.1, the size limit
// of the page) from NCBI through the catalog menu and builds the kallisto index. Not a test: a measurement for
// docs/browser/page-v1-measurements.md, "Page genome path: peak memory". Linux only (/proc).
//   node genome-memory.mjs <absolute path of dist-site> <out.tsv> [engines, default chromium,firefox,webkit] [runs, default 3]
// Each run is a new browser process. The whole-browser Pss is sampled every 250 ms (procTree, sampleTree and Sampler are copied
// from validation/browser-page/measure/measure.mjs, not changed apart from the names of the constants).
import http from 'node:http';
import { readFileSync, readdirSync, statSync, existsSync, writeFileSync } from 'node:fs';
import { join, resolve, extname } from 'node:path';
import { chromium, firefox, webkit } from 'playwright';

const root = resolve(process.argv[2] ?? '');
const outFile = process.argv[3];
if (!process.argv[2] || !outFile || !existsSync(join(root, 'analyze/index.html'))) {
  console.error('usage: node genome-memory.mjs <absolute path of dist-site> <out.tsv> [engines] [runs]');
  process.exit(2);
}
const engines = (process.argv[4] ?? 'chromium,firefox,webkit').split(',').filter(Boolean);
const RUNS = Number(process.argv[5] ?? 3);
const ACCESSION = 'GCF_020509005.1';
const SAMPLE_MS = 250;
const GUARD_KB = 10 * 1024 * 1024; // stop a run above 10 GiB tree Pss (as measure.mjs)
const TIMEOUT_MS = 20 * 60 * 1000;

// ---------------------------------------------------------------- copied from measure.mjs (process tree sampling)
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
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json' };
const server = http.createServer((req, res) => {
  const path = decodeURIComponent(new URL(req.url, 'http://x').pathname);
  let file = resolve(join(root, path));
  if (!file.startsWith(root)) { res.writeHead(403).end(); return; }
  if (existsSync(file) && statSync(file).isDirectory()) file = join(file, 'index.html');
  if (!existsSync(file)) { res.writeHead(404).end('not found'); return; }
  res.writeHead(200, { 'content-type': TYPES[extname(file)] ?? 'application/octet-stream' }).end(readFileSync(file));
});
await new Promise((r) => server.listen(0, '127.0.0.1', r));
const base = `http://127.0.0.1:${server.address().port}`;

const COLS = ['engine', 'browser_version', 'accession', 'run', 'outcome', 'seconds_to_ready', 'index_time', 'genes', 'transcripts',
  'load_tree_pss_kb', 'peak_tree_pss_kb', 'peak_phase', 'peak_tree_rss_sum_kb', 'samples', 'detail'];
const rows = [];
const clean = (v) => String(v ?? '').replace(/[\t\r\n]+/g, ' ').slice(0, 500);

for (const engine of engines) {
  for (let run = 1; run <= RUNS; run++) {
    const row = { engine, accession: ACCESSION, run };
    let browser;
    let sampler;
    try {
      browser = await { chromium, firefox, webkit }[engine].launch();
      row.browser_version = browser.version();
      const page = await browser.newPage();
      await page.goto(`${base}/analyze/`);
      await page.waitForSelector('#step-genome[data-ready="menu"]', { timeout: 60_000 });
      await page.waitForTimeout(500);
      row.load_tree_pss_kb = sampleTree().pss;
      let phase = 'choose';
      let peakPhase = '';
      let guard = null;
      sampler = new Sampler((pss) => { guard = `mem_guard: tree Pss ${pss} kB`; killTree(); });
      const sample = sampler.sample.bind(sampler);
      sampler.sample = () => { const before = sampler.peak?.pss ?? 0; const s = sample(); if (s.pss > before) peakPhase = phase; return s; };
      await page.fill('#genome-filter', ACCESSION);
      await page.selectOption('#genome-menu', ACCESSION);
      const started = Date.now();
      await page.click('#genome-load');
      for (;;) {
        if (guard) throw new Error(guard);
        if (Date.now() - started > TIMEOUT_MS) throw new Error(`timeout after ${TIMEOUT_MS} ms`);
        phase = await page.getAttribute('#step-genome', 'data-phase');
        if (phase === 'ready') break;
        if (!(await page.isHidden('#genome-error'))) throw new Error(`page error message: ${await page.textContent('#genome-error')}`);
        await page.waitForTimeout(250);
      }
      row.seconds_to_ready = ((Date.now() - started) / 1000).toFixed(2);
      row.outcome = 'finished';
      row.index_time = await page.textContent('#ready-time');
      row.genes = await page.textContent('#ready-genes');
      row.transcripts = await page.textContent('#ready-transcripts');
      row.peak_phase = peakPhase;
    } catch (error) {
      row.outcome = /mem_guard/.test(error.message) ? 'mem_guard' : /timeout/.test(error.message) ? 'timeout' : 'error';
      row.detail = error.message.split('\n')[0];
    } finally {
      if (sampler) {
        sampler.stop();
        Object.assign(row, { peak_tree_pss_kb: sampler.peak.pss, peak_tree_rss_sum_kb: sampler.peak.rss, samples: sampler.peak.samples });
      }
      await Promise.race([browser?.close().catch(() => {}), new Promise((done) => setTimeout(done, 15000))]);
      killTree();
    }
    rows.push(row);
    console.log(COLS.map((c) => `${c}=${clean(row[c])}`).join(' '));
  }
}
server.close();
writeFileSync(outFile, [COLS.join('\t'), ...rows.map((r) => COLS.map((c) => clean(r[c])).join('\t'))].join('\n') + '\n');
process.exit(0);
