#!/usr/bin/env node
/**
 * Browser spike 0, Task 3: run the WASM tools on the spike 0 fixture, in Node or in a Playwright browser.
 * Runs inside the stajichlab/genomics-wasm fork (branch qrnamap-spike0), from the repository root.
 * Templates: scripts/validate.mjs (HISAT2), scripts/validate-quant.mjs (kallisto),
 * scripts/validate-rna-analysis.mjs (featureCounts), scripts/validate-deseq2.mjs (DESeq2) at the pin.
 *
 * Usage:
 *   node spike0/parity.mjs --engine node|chromium|firefox|webkit --tools kallisto,hisat2,featurecounts
 *        --data DIR --out DIR [--root DIR] [--deseq2-root DIR] [--rscript FILE]
 *
 * --root         repository root with runtime/ and dist/<tool>/ (default: current directory)
 * --data         staged inputs (plain files, see STAGED below); made by the workflow
 * --out          output directory; one sub-directory per tool and one <tool>.steps.tsv per tool
 * --deseq2-root  unpacked deseq2-webr package (has runtime/deseq2.mjs and dist/deseq2/)
 * --rscript      deseq2_ref.R (same script as the native reference)
 *
 * Each tool reads the same bytes as its native reference (validation/browser-spike0/native_refs.sh):
 *   kallisto      yeast.transcripts.fa (zstd -d of the fixture); reads as .fastq.gz (native read .gz)
 *   HISAT2        yeast.fa (gunzip of the NCBI .fna.gz); reads as plain FASTQ (native read zcat output)
 *   featureCounts yeast.genes.saf; yeast.sam = the native HISAT2 SAM (so only featureCounts is tested)
 *   DESeq2        amb_gene_counts.tsv
 * Tool arguments equal the native arguments, apart from file paths.
 * A failed step is a result: it is written to the steps TSV and the next steps of that tool are skipped.
 * Exit status: 0 when every planned step was attempted or skipped and recorded; 2 on a harness error.
 */
import { readFile, writeFile, mkdir } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { createServer } from 'node:http';
import { dirname, resolve, join, relative, sep } from 'node:path';
import { pathToFileURL } from 'node:url';

const options = { engine: 'node', tools: 'kallisto,hisat2,featurecounts', data: '', out: '', root: process.cwd(), 'deseq2-root': '', rscript: '' };
for (let i = 2; i < process.argv.length; i += 2) {
  const key = process.argv[i].replace(/^--/, '');
  if (!(key in options) || process.argv[i + 1] === undefined) throw new Error(`unknown or empty option ${process.argv[i]}`);
  options[key] = process.argv[i + 1];
}
const ENGINES = ['node', 'chromium', 'firefox', 'webkit'];
if (!ENGINES.includes(options.engine)) throw new Error(`engine must be one of ${ENGINES}`);
if (!options.data || !options.out) throw new Error('--data and --out are required');
const root = resolve(options.root), dataDir = resolve(options.data), outDir = resolve(options.out);
const tools = options.tools.split(',').filter(Boolean);
const timeoutMs = 60 * 60 * 1000; // per command; native times are below 30 s (native_times.tsv)
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const log = message => console.log(`[${new Date().toISOString()}] ${options.engine} ${message}`);
await mkdir(outDir, { recursive: true });

// One plan per tool. Input src 'data:NAME' is a staged input; 'out:PATH' is an output of an earlier step.
const HT2 = Array.from({ length: 8 }, (_, i) => `yeast.hs.${i + 1}.ht2`);
const PLANS = {
  kallisto: { waitForExit: false, steps: [
    { step: 'index', args: ['index', '-i', 'yeast.idx', 'yeast.transcripts.fa'],
      inputs: [['yeast.transcripts.fa', 'data:yeast.transcripts.fa']], outputs: ['yeast.idx'] },
    { step: 'quant', args: ['quant', '-i', 'yeast.idx', '-o', 'kq', '-t', '1', '--plaintext', 'yeast_sub_1.fastq.gz', 'yeast_sub_2.fastq.gz'],
      inputs: [['yeast.idx', 'out:yeast.idx'], ['yeast_sub_1.fastq.gz', 'data:yeast_sub_1.fastq.gz'], ['yeast_sub_2.fastq.gz', 'data:yeast_sub_2.fastq.gz']],
      outputs: ['kq/abundance.tsv', 'kq/run_info.json'] },
  ] },
  hisat2: { waitForExit: false, steps: [
    { step: 'build', args: ['build', '-p', '1', 'yeast.fa', 'yeast.hs'], inputs: [['yeast.fa', 'data:yeast.fa']], outputs: HT2 },
    { step: 'align', args: ['align', '-p', '1', '-x', 'yeast.hs', '-1', 'reads_1.fastq', '-2', 'reads_2.fastq', '-S', 'yeast.sam'],
      inputs: [...HT2.map(f => [f, `out:${f}`]), ['reads_1.fastq', 'data:reads_1.fastq'], ['reads_2.fastq', 'data:reads_2.fastq']],
      outputs: ['yeast.sam'] },
  ] },
  featurecounts: { waitForExit: true, steps: [
    { step: 'count', args: ['-T', '1', '-F', 'SAF', '-a', 'yeast.genes.saf', '-p', '--countReadPairs', '-s', '0', '-o', 'fc.tsv', 'yeast.sam'],
      inputs: [['yeast.genes.saf', 'data:yeast.genes.saf'], ['yeast.sam', 'data:native_yeast.sam']], outputs: ['fc.tsv', 'fc.tsv.summary'] },
  ] },
};

const meta = { engine: options.engine, node: process.version, startedAt: new Date().toISOString(), root, tools, builds: {}, inputs: {} };
let browser, server, page;

function srcPath(tool, src) {
  const [kind, name] = [src.slice(0, src.indexOf(':')), src.slice(src.indexOf(':') + 1)];
  return kind === 'data' ? join(dataDir, name) : join(outDir, tool, name);
}

function inside(base, path) { return path === base || path.startsWith(base + sep); }

async function startBrowser() {
  const playwright = await import('playwright'); // resolved from the repository node_modules (npm ci)
  browser = await playwright[options.engine].launch({ headless: true });
  meta.browserVersion = browser.version();
  server = createServer(async (request, response) => {
    // Cross-origin isolation, as in the upstream validation servers (featureCounts needs pthreads).
    response.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
    response.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
    response.setHeader('Cross-Origin-Resource-Policy', 'same-origin');
    try {
      const url = new URL(request.url, 'http://localhost');
      const name = decodeURIComponent(url.pathname);
      if (name === '/') { response.setHeader('Content-Type', 'text/html'); response.end('<!doctype html><meta charset="utf-8"><title>spike0 parity</title>'); return; }
      if (name.startsWith('/out/')) {
        const path = resolve(outDir, '.' + name.slice(4));
        if (!inside(outDir, path)) { response.writeHead(403).end(); return; }
        if (request.method === 'PUT') {
          const chunks = [];
          for await (const chunk of request) chunks.push(chunk);
          await mkdir(dirname(path), { recursive: true });
          await writeFile(path, Buffer.concat(chunks));
          response.writeHead(204).end();
          return;
        }
        response.setHeader('Content-Type', 'application/octet-stream');
        response.end(await readFile(path));
        return;
      }
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
  await new Promise(done => server.listen(0, '127.0.0.1', done));
  page = await browser.newPage();
  page.on('console', message => log(`browser ${message.type()}: ${message.text()}`));
  page.on('pageerror', error => log(`browser pageerror: ${error.message}`));
  page.on('crash', () => log('browser page crashed'));
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  meta.crossOriginIsolated = await page.evaluate(() => crossOriginIsolated);
  log(`browser ${meta.browserVersion}, crossOriginIsolated=${meta.crossOriginIsolated}`);
}

// Node: read inputs from disk, run, write outputs. Returns { exitCode, missingOutputs, files: [paths] }.
async function runNode(tool, plan, step) {
  const { ToolRunner } = await import(pathToFileURL(join(root, 'runtime', 'index.mjs')).href);
  const files = [];
  for (const [path, src] of step.inputs) files.push({ path, data: await readFile(srcPath(tool, src)) });
  const runner = new ToolRunner({ moduleUrl: pathToFileURL(join(root, 'dist', tool, `${tool}.mjs`)), waitForExit: plan.waitForExit });
  const result = await runner.run(step.args, { files, outputs: step.outputs, timeoutMs });
  const dir = join(outDir, tool);
  await mkdir(dir, { recursive: true });
  await writeFile(join(dir, `${step.step}.stdout.txt`), result.stdout);
  await writeFile(join(dir, `${step.step}.stderr.txt`), result.stderr);
  const written = [];
  for (const file of result.files) {
    const rel = file.path.replace(/^\/work\//, '');
    await mkdir(dirname(join(dir, rel)), { recursive: true });
    await writeFile(join(dir, rel), file.data);
    written.push(rel);
  }
  return { exitCode: result.exitCode, missingOutputs: result.missingOutputs, files: written };
}

// Browser: the page fetches inputs from this server, runs ToolRunner in a worker, and PUTs outputs back.
async function runBrowser(tool, plan, step) {
  const inputs = step.inputs.map(([path, src]) => {
    const kind = src.slice(0, src.indexOf(':')), name = src.slice(src.indexOf(':') + 1);
    return { path, url: kind === 'data' ? `/data/${name}` : `/out/${tool}/${name}` };
  });
  const evaluation = page.evaluate(async ({ tool, stepName, args, inputs, outputs, waitForExit, timeoutMs }) => {
    const { ToolRunner } = await import('/runtime/index.mjs');
    const put = async (name, data) => {
      const response = await fetch(`/out/${tool}/${name}`, { method: 'PUT', body: data });
      if (!response.ok) throw new Error(`upload ${name}: HTTP ${response.status}`);
    };
    const files = [];
    for (const input of inputs) {
      const response = await fetch(input.url);
      if (!response.ok) throw new Error(`fetch ${input.url}: HTTP ${response.status}`);
      files.push({ path: input.path, data: new Uint8Array(await response.arrayBuffer()) });
    }
    const runner = new ToolRunner({ moduleUrl: new URL(`/dist/${tool}/${tool}.mjs`, location.href), waitForExit });
    const result = await runner.run(args, { files, outputs, timeoutMs });
    files.length = 0;
    await put(`${stepName}.stdout.txt`, result.stdout);
    await put(`${stepName}.stderr.txt`, result.stderr);
    const written = [];
    for (const file of result.files) {
      const rel = file.path.replace(/^\/work\//, '');
      await put(rel, file.data);
      written.push(rel);
    }
    return { exitCode: result.exitCode, missingOutputs: result.missingOutputs, files: written };
  }, { tool, stepName: step.step, args: step.args, inputs, outputs: step.outputs, waitForExit: plan.waitForExit, timeoutMs });
  let timer;
  try {
    return await Promise.race([evaluation, new Promise((_, reject) => {
      timer = setTimeout(() => reject(new Error(`browser step exceeded ${timeoutMs + 300000} ms`)), timeoutMs + 300000);
    })]);
  } finally { clearTimeout(timer); }
}

async function writeSteps(tool, rows) {
  const header = ['engine', 'tool', 'step', 'status', 'exit_code', 'seconds', 'args', 'detail', 'outputs'];
  const clean = value => String(value ?? '').replace(/[\t\r\n]+/g, ' ');
  const text = [header.join('\t'), ...rows.map(row => header.map(key => clean(row[key])).join('\t'))].join('\n') + '\n';
  await writeFile(join(outDir, `${tool}.steps.tsv`), text);
}

async function runTool(tool) {
  const plan = PLANS[tool];
  const rows = [];
  let failed = false;
  for (const step of plan.steps) {
    const row = { engine: options.engine, tool, step: step.step, args: step.args.join(' ') };
    rows.push(row);
    if (failed) { Object.assign(row, { status: 'skipped', detail: 'an earlier step of this tool did not succeed' }); continue; }
    log(`${tool} ${step.step} starting: ${row.args}`);
    const started = Date.now();
    try {
      const result = options.engine === 'node' ? await runNode(tool, plan, step) : await runBrowser(tool, plan, step);
      row.seconds = ((Date.now() - started) / 1000).toFixed(2);
      row.exit_code = result.exitCode;
      const outputs = [];
      for (const rel of result.files) outputs.push(`${rel}=${hash(await readFile(join(outDir, tool, rel)))}`);
      row.outputs = outputs.join(';');
      if (result.exitCode !== 0 || result.missingOutputs.length) {
        let tail = '';
        try { tail = (await readFile(join(outDir, tool, `${step.step}.stderr.txt`), 'utf8')).slice(-600); } catch {}
        Object.assign(row, { status: 'fail', detail: `exit ${result.exitCode}; missing ${result.missingOutputs.join(',') || 'none'}; stderr tail: ${tail}` });
        failed = true;
      } else row.status = 'ok';
    } catch (error) {
      row.seconds = ((Date.now() - started) / 1000).toFixed(2);
      Object.assign(row, { status: 'error', detail: String(error?.stack || error).slice(0, 1500) });
      failed = true;
      if (page && page.isClosed()) { log('page closed; opening a new page'); page = await browser.newPage(); await page.goto(`http://127.0.0.1:${server.address().port}/`); }
    }
    log(`${tool} ${step.step} ${row.status} exit=${row.exit_code ?? ''} after ${row.seconds} s ${row.status === 'ok' ? '' : row.detail}`);
    await writeSteps(tool, rows);
  }
  await writeSteps(tool, rows);
}

// DESeq2 in Node only: the native deseq2_ref.R, unchanged, with commandArgs() and the output directory set by a prelude.
async function runDESeq2() {
  const rows = [];
  const row = { engine: options.engine, tool: 'deseq2', step: 'wald', args: 'deseq2_ref.R input/amb_gene_counts.tsv output/amb' };
  rows.push(row);
  const started = Date.now();
  let session;
  try {
    if (options.engine !== 'node') throw new Error('DESeq2 is run in Node only');
    const pkg = resolve(options['deseq2-root']);
    const runtimeFile = join(pkg, 'runtime', 'deseq2.mjs');
    meta.builds.deseq2 = { runtime: hash(await readFile(runtimeFile)), manifest: hash(await readFile(join(pkg, 'dist', 'deseq2', 'manifest.json'))) };
    const manifest = JSON.parse(await readFile(join(pkg, 'dist', 'deseq2', 'manifest.json'), 'utf8'));
    meta.builds.deseq2.versions = { deseq2: manifest.version, r: manifest.rVersion, webr: manifest.webrVersion };
    const script = await readFile(resolve(options.rscript), 'utf8');
    meta.builds.deseq2.script = hash(script);
    const prelude = 'commandArgs <- function(trailingOnly = FALSE) c("input/amb_gene_counts.tsv", "output/amb")\n'
      + 'dir.create("output", showWarnings = FALSE)\n';
    const counts = await readFile(join(dataDir, 'amb_gene_counts.tsv'));
    meta.inputs['amb_gene_counts.tsv'] = { bytes: counts.length, sha256: hash(counts) };
    const { createDESeq2 } = await import(pathToFileURL(runtimeFile).href);
    session = await createDESeq2({ onProgress: message => log(`deseq2 ${message}`) });
    const result = await session.run(prelude + script, { files: { 'input/amb_gene_counts.tsv': counts } });
    const dir = join(outDir, 'deseq2');
    await mkdir(dir, { recursive: true });
    await writeFile(join(dir, 'wald.stdout.txt'), result.stdout || '');
    await writeFile(join(dir, 'wald.stderr.txt'), result.stderr || '');
    const outputs = [];
    for (const [path, bytes] of Object.entries(result.files)) {
      if (!path.startsWith('output/')) continue;
      const name = path.slice('output/'.length);
      if (!name || name.includes('/')) continue;
      await writeFile(join(dir, name), new Uint8Array(bytes));
      outputs.push(`${name}=${hash(new Uint8Array(bytes))}`);
    }
    row.outputs = outputs.join(';');
    row.exit_code = 0;
    row.status = outputs.some(o => o.startsWith('amb.results.csv=')) ? 'ok' : 'fail';
    if (row.status !== 'ok') row.detail = `no amb.results.csv; stderr: ${String(result.stderr).slice(-600)}`;
  } catch (error) {
    Object.assign(row, { status: 'error', detail: String(error?.stack || error).slice(0, 1500) + (error?.stderr ? ` R stderr: ${String(error.stderr).slice(-600)}` : '') });
  } finally {
    try { session?.close(); } catch {}
    row.seconds = ((Date.now() - started) / 1000).toFixed(2);
    log(`deseq2 wald ${row.status} after ${row.seconds} s ${row.detail ?? ''}`);
    await writeSteps('deseq2', rows);
  }
}

let harnessError = false;
try {
  for (const tool of tools) {
    if (tool === 'deseq2') continue;
    if (!PLANS[tool]) throw new Error(`unknown tool ${tool}`);
    const assets = {};
    for (const name of [`${tool}.mjs`, `${tool}.wasm`, `${tool}.worker.mjs`]) {
      try { assets[name] = hash(await readFile(join(root, 'dist', tool, name))); } catch {}
    }
    meta.builds[tool] = assets;
    for (const step of PLANS[tool].steps) for (const [, src] of step.inputs) {
      if (!src.startsWith('data:')) continue;
      const name = src.slice(5);
      if (!meta.inputs[name]) { const bytes = await readFile(join(dataDir, name)); meta.inputs[name] = { bytes: bytes.length, sha256: hash(bytes) }; }
    }
  }
  meta.runtime = {};
  for (const name of ['index.mjs', 'worker.mjs', 'execute.mjs']) meta.runtime[name] = hash(await readFile(join(root, 'runtime', name)));
  if (options.engine !== 'node' && tools.some(t => t !== 'deseq2')) await startBrowser();
  for (const tool of tools) {
    if (tool === 'deseq2') await runDESeq2();
    else await runTool(tool);
  }
} catch (error) {
  harnessError = true;
  meta.harnessError = String(error?.stack || error);
  console.error(error?.stack || error);
} finally {
  meta.finishedAt = new Date().toISOString();
  await writeFile(join(outDir, `meta-${tools.join('-')}.json`), JSON.stringify(meta, null, 2) + '\n');
  if (browser) await browser.close().catch(() => {});
  if (server) await new Promise(done => server.close(done));
}
process.exit(harnessError ? 2 : 0);
