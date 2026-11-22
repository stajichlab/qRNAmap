// Browser smoke test: node smoke.mjs <absolute path of dist-site> [engines, comma list; default chromium,firefox,webkit]
// Serves dist-site with a small static server, loads / and /analyze/ in each engine, and fails on a console error,
// a failed request, a page error, or a broken "Analyze reads" link (it clicks the button and checks the next page).
import http from 'node:http';
import { readFileSync, statSync, existsSync } from 'node:fs';
import { join, resolve, extname } from 'node:path';
import { chromium, firefox, webkit } from 'playwright';

const root = resolve(process.argv[2] ?? '');
if (!process.argv[2] || !existsSync(root)) { console.error('usage: node smoke.mjs <absolute path of dist-site> [engines]'); process.exit(2); }
const engines = (process.argv[3] ?? 'chromium,firefox,webkit').split(',').filter(Boolean);
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json', '.md': 'text/plain', '.txt': 'text/plain' };

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

const launchers = { chromium, firefox, webkit };
const failures = [];
for (const name of engines) {
  let browser;
  try {
    browser = await launchers[name].launch();
    const page = await browser.newPage({ viewport: { width: 390, height: 800 } });
    const errors = [];
    page.on('console', (m) => { if (m.type() === 'error') errors.push(`console error: ${m.text()}`); });
    page.on('pageerror', (e) => errors.push(`page error: ${e.message}`));
    page.on('requestfailed', (r) => errors.push(`request failed: ${r.url()} ${r.failure()?.errorText}`));
    page.on('response', (r) => { if (r.status() >= 400) errors.push(`HTTP ${r.status()}: ${r.url()}`); });

    await page.goto(`${base}/`, { waitUntil: 'load' });
    if (!(await page.textContent('body')).includes('Status: only the genome step of the analysis tool works')) errors.push('landing status line missing');
    const href = await page.getAttribute('a.button', 'href');
    if (href !== 'analyze/') errors.push(`button href is ${href}`);
    // no horizontal scroll at phone width
    const wide = await page.evaluate(() => document.documentElement.scrollWidth > document.documentElement.clientWidth);
    if (wide) errors.push('/ scrolls horizontally at 390 px');
    await Promise.all([page.waitForURL(`${base}/analyze/`), page.click('a.button')]);
    await page.waitForLoadState('load');
    await page.waitForSelector('#step-genome[data-ready="menu"]', { timeout: 30000 }); // the genome menu has loaded (wasm catalog)
    // the three steps and their enable order
    const states = await page.$$eval('.step', (els) => els.map((e) => `${e.dataset.step}:${e.dataset.state}:${e.querySelector('fieldset').disabled}`));
    if (states.join(',') !== 'genome:active:false,samples:locked:true,results:locked:true') errors.push(`step states: ${states.join(',')}`);
    await page.evaluate(() => window.qrnamapShell.completeStep('genome'));
    const after = await page.$$eval('.step', (els) => els.map((e) => `${e.dataset.state}:${e.querySelector('fieldset').disabled}`));
    if (after.join(',') !== 'done:false,active:false,locked:true') errors.push(`step states after genome done: ${after.join(',')}`);
    if (!(await page.textContent('#construction-notice')).includes('Under construction: genome step works; samples and results are not built yet.')) errors.push('construction notice missing');
    if (await page.$('#preview-label')) errors.push('preview label present before a tool state exists');
    if (await page.evaluate(() => document.documentElement.scrollWidth > document.documentElement.clientWidth)) errors.push('/analyze/ scrolls horizontally at 390 px');
    console.log(`${name}: ${errors.length ? 'FAIL' : 'ok'} (${await browser.version()})`);
    for (const e of errors) { console.log(`  ${e}`); failures.push(`${name}: ${e}`); }
  } catch (e) {
    console.log(`${name}: FAIL ${e.message.split('\n')[0]}`);
    failures.push(`${name}: ${e.message.split('\n')[0]}`);
  } finally {
    await browser?.close();
  }
}
server.close();
process.exit(failures.length ? 1 : 0);
