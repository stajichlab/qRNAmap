// Tests for the link checker, run against temp copies of a small fixture dist: node --test validation/browser-page/smoke/
import test from 'node:test';
import assert from 'node:assert/strict';
import { mkdtempSync, mkdirSync, writeFileSync, rmSync, readFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { checkDist } from './links-lib.mjs';

const FILES = {
  'index.html': '<!doctype html><link rel="stylesheet" href="style.css?v=1"><a class="button" href="analyze/">Go</a>',
  'analyze/index.html': '<!doctype html><script type="module" src="app.js?v=1"></script><link rel="stylesheet" href="../style.css">',
  'analyze/app.js': "import { a } from '../lib.js?v=1';\nexport { a } from '../lib.js';\nimport '../lib.js';\nconst m = () => import('../lib.js');\n",
  'lib.js': 'export const a = 1;\n',
  'style.css': '@import "base.css";\nbody { background: url("bg.png"); }\n.x { background: url(data:image/png;base64,AAAA); }\n',
  'base.css': 'h1 { color: red; }\n',
  'bg.png': 'x',
};
function fixture(mutate) {
  const dir = mkdtempSync(join(tmpdir(), 'dist-fixture-'));
  for (const [name, text] of Object.entries(FILES)) {
    mkdirSync(join(dir, name, '..'), { recursive: true });
    writeFileSync(join(dir, name), text);
  }
  mutate?.(dir);
  return dir;
}
const edit = (dir, name, f) => writeFileSync(join(dir, name), f(readFileSync(join(dir, name), 'utf8')));
function problemsAfter(mutate) {
  const dir = fixture(mutate);
  try { return checkDist(dir); } finally { rmSync(dir, { recursive: true, force: true }); }
}

test('the fixture passes', () => assert.deepEqual(problemsAfter(), []));
test('a missing static import target fails', () => {
  const p = problemsAfter((d) => edit(d, 'analyze/app.js', (t) => t.replace("from '../lib.js?v=1'", "from '../gone.js?v=1'")));
  assert.ok(p.some((x) => x.includes('gone.js')), p.join('\n'));
});
test('a missing dynamic import target fails', () => {
  const p = problemsAfter((d) => edit(d, 'analyze/app.js', (t) => t.replace("import('../lib.js')", "import('../gone2.js')")));
  assert.ok(p.some((x) => x.includes('gone2.js')), p.join('\n'));
});
test('a missing export-from target fails', () => {
  const p = problemsAfter((d) => edit(d, 'analyze/app.js', (t) => t.replace("export { a } from '../lib.js'", "export { a } from '../gone3.js'")));
  assert.ok(p.some((x) => x.includes('gone3.js')), p.join('\n'));
});
test('a missing bare import target fails', () => {
  const p = problemsAfter((d) => edit(d, 'analyze/app.js', (t) => t.replace("import '../lib.js';", "import '../gone4.js';")));
  assert.ok(p.some((x) => x.includes('gone4.js')), p.join('\n'));
});
test('a missing CSS url target fails', () => {
  const p = problemsAfter((d) => edit(d, 'style.css', (t) => t.replace('bg.png', 'gone.png')));
  assert.ok(p.some((x) => x.includes('gone.png')), p.join('\n'));
});
test('a missing CSS @import target fails', () => {
  const p = problemsAfter((d) => edit(d, 'style.css', (t) => t.replace('base.css', 'gone.css')));
  assert.ok(p.some((x) => x.includes('gone.css')), p.join('\n'));
});
test('a missing <script src> fails', () => {
  const p = problemsAfter((d) => edit(d, 'analyze/index.html', (t) => t.replace('app.js', 'gone.js')));
  assert.ok(p.some((x) => x.includes('gone.js')), p.join('\n'));
});
test('a broken button href fails', () => {
  const p = problemsAfter((d) => edit(d, 'index.html', (t) => t.replace('href="analyze/"', 'href="analyse/"')));
  assert.ok(p.some((x) => x.includes('button')), p.join('\n'));
});
test('a leftover __BUILD__ token fails', () => {
  const p = problemsAfter((d) => edit(d, 'lib.js', (t) => t + "// ?v=__BUILD__\n"));
  assert.ok(p.some((x) => x.includes('__BUILD__')), p.join('\n'));
});
