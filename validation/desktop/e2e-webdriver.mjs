// End-to-end check of the desktop app (Linux): the built app is started by tauri-driver (WebKitWebDriver) and driven
// through its own page. No npm packages: plain W3C WebDriver over HTTP.
//
//   APP=/usr/bin/qrnamap-desktop QRNAMAP_KALLISTO=/path/to/kallisto \
//     node validation/desktop/e2e-webdriver.mjs <reads_1.fastq.gz> <reads_2.fastq.gz> [accession]
//
// Needs: tauri-driver (cargo install tauri-driver --locked), WebKitWebDriver (Debian/Ubuntu: webkit2gtk-driver), a display
// (xvfb-run on a headless machine), network access to NCBI and the S3 index mirror. Reads: yeast SRR3396381 from release
// testsets-v1 fit the default accession (S. cerevisiae S288C, GCF_000146045.2).
//
// Checks, in order:
//   1. the page sees the Tauri bridge (html.desktop) and check_native_capabilities finds kallisto;
//   2. the genome step of the page loads the accession and leaves a native index that the native kallisto loads (on Linux
//      the bioconda kallisto cannot read the page's kallisto.wasm index, so this is the "build a native index" path);
//   3. run_native_quant counts the reads with that index: every pair processed, more than half pseudoaligned.
// Exit 0 when all pass. Each step prints one line.
import { spawn } from 'node:child_process';
import { resolve } from 'node:path';
import { existsSync } from 'node:fs';

const [r1, r2, accession = 'GCF_000146045.2'] = process.argv.slice(2);
if (!r1 || !r2) {
  console.error('usage: node e2e-webdriver.mjs <reads_1.fastq.gz> <reads_2.fastq.gz> [accession]');
  process.exit(2);
}
for (const f of [r1, r2]) if (!existsSync(f)) { console.error(`no file ${f}`); process.exit(2); }
const APP = process.env.APP || '/usr/bin/qrnamap-desktop';
const PORT = Number(process.env.WEBDRIVER_PORT || 4444);
const GENOME_MS = Number(process.env.GENOME_TIMEOUT_MS || 20 * 60_000);
const WD = `http://127.0.0.1:${PORT}`;

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const t0 = Date.now();
const log = (msg) => console.log(`[${((Date.now() - t0) / 1000).toFixed(1).padStart(6)} s] ${msg}`);

async function wd(method, path, body) {
  const res = await fetch(`${WD}${path}`, { method, headers: { 'content-type': 'application/json' }, body: body ? JSON.stringify(body) : undefined });
  const json = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(`${method} ${path}: HTTP ${res.status} ${JSON.stringify(json.value ?? json).slice(0, 500)}`);
  return json.value;
}

const driver = spawn(process.env.TAURI_DRIVER || 'tauri-driver', ['--port', String(PORT)], { stdio: ['ignore', 'inherit', 'inherit'] });
let session = null;
let failed = false;
const fail = (msg) => { failed = true; log(`FAIL ${msg}`); };

try {
  for (let i = 0; ; i++) {
    try { await fetch(`${WD}/status`); break; } catch { if (i > 50) throw new Error('tauri-driver did not start'); await sleep(200); }
  }
  session = (await wd('POST', '/session', { capabilities: { alwaysMatch: { 'tauri:options': { application: APP } } } })).sessionId;
  const S = `/session/${session}`;
  const exec = (script, args = []) => wd('POST', `${S}/execute/sync`, { script, args });
  // An async script resolves through its last argument; errors come back as { error }.
  const execAsync = (body, args = []) => wd('POST', `${S}/execute/async`, {
    script: `const done = arguments[arguments.length - 1]; (async () => { ${body} })().then((v) => done({ v }), (e) => done({ error: String(e?.message || e) }));`,
    args,
  });
  const invoke = async (cmd, args) => {
    const r = await execAsync('return await window.__TAURI_INTERNALS__.invoke(arguments[0], arguments[1]);', [cmd, args]);
    if (r.error) throw new Error(`${cmd}: ${r.error}`);
    return r.v;
  };
  await wd('POST', `${S}/timeouts`, { script: 30 * 60_000 });

  const start = await wd('GET', `${S}/url`);
  await wd('POST', `${S}/url`, { url: new URL('analyze/', start).href });
  for (let i = 0; !(await exec('return Boolean(window.qrnamapApp && window.qrnamapApp.loadCatalogAccession)')); i++) {
    if (i > 150) throw new Error('the /analyze/ page did not start (no window.qrnamapApp)');
    await sleep(200);
  }
  log(`app started, page ${await wd('GET', `${S}/url`)}`);

  // 1. Tauri bridge and native kallisto
  const desktop = await exec("return document.documentElement.classList.contains('desktop')");
  const caps = await invoke('check_native_capabilities', {});
  log(`1. html.desktop=${desktop}; native: ${JSON.stringify(caps)}`);
  if (!desktop) fail('the page did not mark itself as the desktop app');
  if (!caps?.available) throw new Error('check_native_capabilities: kallisto not available (set QRNAMAP_KALLISTO?)');

  // 2. Genome step through the page, then a native index
  await exec('window.qrnamapApp.loadCatalogAccession(arguments[0]).catch((e) => { window.__e2eGenomeError = String(e?.message || e); });', [accession]);
  let st;
  for (;;) {
    st = await exec(`const s = window.qrnamapApp.state; return { phase: s.phase, ready: s.indexReady, native: s.nativeIndexPath,
      nativeErr: s.nativeIndexError, note: s.nativeIndexNote, genes: s.nGenes, tx: s.nTranscripts, err: window.__e2eGenomeError || null,
      shown: (document.getElementById('genome-error') || {}).textContent || '' }`);
    if (st.err) throw new Error(`genome step: ${st.err}`);
    if (st.phase === 'ready' && st.ready) break;
    if (Date.now() - t0 > GENOME_MS) throw new Error(`genome step did not finish in ${GENOME_MS / 1000} s: ${JSON.stringify(st)}`);
    await sleep(2000);
  }
  log(`2. genome ${accession}: ${st.genes} genes, ${st.tx} transcripts; native index ${st.native || 'none'}${st.note ? ` (${st.note})` : ''}`);
  if (!st.native) throw new Error(`no native index: ${st.nativeErr || 'unknown reason'}`);
  const check = await invoke('check_native_index', { indexPath: st.native });
  if (!check.ok) fail(`the native kallisto does not load ${st.native}: ${check.detail}`);

  // 3. Native quantification of the reads
  const tx2gene = await exec('return window.qrnamapApp.state.tx2gene');
  const q = await invoke('run_native_quant', { req: {
    name: 'e2e', fastq1: [resolve(r1)], fastq2: [resolve(r2)], index_path: st.native, tx2gene,
    strandedness: 'unstranded', fragment_mean: null, fragment_sd: null, threads: null, out_dir: null,
  } });
  const genesCounted = Object.values(q.gene_counts || {}).filter((n) => n > 0).length;
  const rate = q.n_processed ? q.n_pseudoaligned / q.n_processed : 0;
  log(`3. native quant: ${q.n_processed} processed, ${q.n_pseudoaligned} pseudoaligned (${(100 * rate).toFixed(1)} %), ${genesCounted} genes with reads, whole file ${q.whole_file}`);
  if (!(q.n_processed > 0)) fail('no reads processed');
  if (!(rate > 0.5)) fail('half or fewer of the reads pseudoaligned');
  if (genesCounted === 0) fail('no gene has a count');
} catch (e) {
  fail(e.message || String(e));
} finally {
  if (session) await wd('DELETE', `/session/${session}`).catch(() => {});
  driver.kill();
}
log(failed ? 'RESULT: FAIL' : 'RESULT: PASS');
process.exit(failed ? 1 : 0);
