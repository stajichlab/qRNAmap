// Peak memory and time of the real quantification path of /analyze/: the page downloads a catalog genome from NCBI through the
// menu and builds the kallisto index; then ONE sample is added through the file chooser (local files, not fetched by the page)
// and counted with the default preview size (N = 200,000 pairs). Not a test: a measurement for docs/browser/page-v1-measurements.md,
// "Page quantification path: peak memory". Linux only (/proc).
//   node quant-memory.mjs <absolute path of dist-site> <out.tsv> <engines> <runs> <case>...
//   case = <label>:<accession>:<mate 1 path>:<mate 2 path>   (absolute paths; the label names the reads in the TSV)
// Each run is a new browser process. The whole-browser Pss is sampled every 250 ms (procTree, sampleTree and Sampler are copied
// from validation/browser-page/measure/measure.mjs through genome-memory.mjs, not changed). The peak is kept per phase:
// genome (download, transcripts, index), prep (reading and cutting the reads), quant (kallisto quant and the gene counts).
import http from 'node:http';
import { readFileSync, readdirSync, statSync, existsSync, writeFileSync } from 'node:fs';
import { join, resolve, extname } from 'node:path';
import { chromium, firefox, webkit } from 'playwright';

const root = resolve(process.argv[2] ?? '');
const outFile = process.argv[3];
const engines = (process.argv[4] ?? '').split(',').filter(Boolean);
const RUNS = Number(process.argv[5]);
const cases = process.argv.slice(6).map((c) => { const [label, accession, mate1, mate2] = c.split(':'); return { label, accession, mate1, mate2 }; });
if (!process.argv[2] || !outFile || !existsSync(join(root, 'analyze/index.html')) || !engines.length || !(RUNS > 0) || !cases.length ||
    cases.some((c) => !c.accession || !existsSync(c.mate1) || !existsSync(c.mate2))) {
  console.error('usage: node quant-memory.mjs <absolute path of dist-site> <out.tsv> <engines> <runs> <label>:<accession>:<mate1>:<mate2>...');
  process.exit(2);
}
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

const PHASES = ['genome', 'prep', 'quant'];
const COLS = ['engine', 'browser_version', 'case', 'accession', 'mate1_bytes', 'mate2_bytes', 'run', 'outcome', 'seconds_genome', 'seconds_sample',
  'pairs', 'processed', 'pseudoaligned', 'reruns', 'load_tree_pss_kb', ...PHASES.map((p) => `${p}_tree_pss_peak_kb`), 'peak_tree_rss_sum_kb', 'samples', 'detail'];
const rows = [];
const clean = (v) => String(v ?? '').replace(/[\t\r\n]+/g, ' ').slice(0, 500);

for (const c of cases) {
  for (const engine of engines) {
    for (let run = 1; run <= RUNS; run++) {
      const row = { engine, case: c.label, accession: c.accession, run, mate1_bytes: statSync(c.mate1).size, mate2_bytes: statSync(c.mate2).size };
      let browser;
      let sampler;
      const peaks = Object.fromEntries(PHASES.map((p) => [p, 0]));
      try {
        browser = await { chromium, firefox, webkit }[engine].launch();
        row.browser_version = browser.version();
        const page = await browser.newPage();
        await page.goto(`${base}/analyze/`);
        await page.waitForSelector('#step-samples[data-ready="table"]', { timeout: 60_000 });
        await page.waitForTimeout(500);
        row.load_tree_pss_kb = sampleTree().pss;
        let phase = 'genome';
        let guard = null;
        sampler = new Sampler((pss) => { guard = `mem_guard: tree Pss ${pss} kB`; killTree(); });
        const sample = sampler.sample.bind(sampler);
        sampler.sample = () => { const s = sample(); if (peaks[phase] !== undefined) peaks[phase] = Math.max(peaks[phase], s.pss); return s; };
        const wait = async (done, what) => {
          const started = Date.now();
          for (;;) {
            if (guard) throw new Error(guard);
            if (Date.now() - started > TIMEOUT_MS) throw new Error(`timeout after ${TIMEOUT_MS} ms (${what})`);
            if (await done()) return;
            await page.waitForTimeout(100);
          }
        };
        await page.fill('#genome-filter', c.accession);
        await page.selectOption('#genome-menu', c.accession);
        let started = Date.now();
        await page.click('#genome-load');
        await wait(async () => {
          if (!(await page.isHidden('#genome-error'))) throw new Error(`genome step error: ${await page.textContent('#genome-error')}`);
          return (await page.getAttribute('#step-genome', 'data-phase')) === 'ready';
        }, 'genome');
        row.seconds_genome = ((Date.now() - started) / 1000).toFixed(2);
        const r = page.locator('#sample-table tbody tr').first();
        await r.locator('.sample-name').fill('S1');
        await r.locator('.cov').first().fill('x');
        await r.locator('td[data-mate="1"] input').setInputFiles(c.mate1);
        await page.locator('#sample-table tbody tr').first().locator('td[data-mate="2"] input').setInputFiles(c.mate2);
        phase = 'prep';
        started = Date.now();
        await page.click('#run-samples');
        await wait(async () => {
          if (!(await page.isHidden('#samples-errors'))) throw new Error(`samples step error: ${await page.textContent('#samples-errors')}`);
          const st = await page.evaluate(() => [window.qrnamapApp.state.samplesPhase, document.getElementById('samples-status').textContent]);
          if (/Running kallisto/.test(st[1])) phase = 'quant';
          else if (/Reading/.test(st[1])) phase = 'prep';
          return ['done', 'cancelled', 'lost'].includes(st[0]);
        }, 'sample');
        row.seconds_sample = ((Date.now() - started) / 1000).toFixed(2);
        const s = await page.evaluate(() => window.qrnamapApp.state.samples[0]);
        Object.assign(row, { pairs: s.nPairs, processed: s.processed, pseudoaligned: s.pseudoaligned, reruns: s.reruns });
        if (s.status !== 'done') throw new Error(`sample ${s.status}: ${s.error}`);
        row.outcome = 'finished';
      } catch (error) {
        row.outcome = /mem_guard/.test(error.message) ? 'mem_guard' : /timeout/.test(error.message) ? 'timeout' : 'error';
        row.detail = error.message.split('\n')[0];
      } finally {
        if (sampler) {
          sampler.stop();
          Object.assign(row, { peak_tree_rss_sum_kb: sampler.peak.rss, samples: sampler.peak.samples });
          for (const p of PHASES) row[`${p}_tree_pss_peak_kb`] = peaks[p];
        }
        await Promise.race([browser?.close().catch(() => {}), new Promise((done) => setTimeout(done, 15000))]);
        killTree();
      }
      rows.push(row);
      console.log(COLS.map((col) => `${col}=${clean(row[col])}`).join(' '));
    }
  }
}
server.close();
writeFileSync(outFile, [COLS.join('\t'), ...rows.map((r) => COLS.map((col) => clean(r[col])).join('\t'))].join('\n') + '\n');
process.exit(0);
