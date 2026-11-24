// Link checks for the built site. checkDist(root) returns a list of problem texts (empty: all good).
// Checks: local href/src in HTML; ES module specifiers in .js/.mjs (import ... from, import '...', import('...'), export ... from);
// CSS url(...) and @import; a leftover __BUILD__ token; the "Analyze reads" button on the landing page. Needs no packages.
import { readFileSync, readdirSync, statSync, existsSync } from 'node:fs';
import { join, resolve, dirname, extname, relative } from 'node:path';

const TEXT = new Set(['.html', '.js', '.mjs', '.css', '.json', '.md', '.txt']);
const walk = (dir) => readdirSync(dir).flatMap((n) => {
  const p = join(dir, n);
  return statSync(p).isDirectory() ? walk(p) : [p];
});
const skip = (v) => v === '' || /^(https?:|mailto:|data:|blob:|#)/i.test(v);

/** Does the local reference `value`, made in file `from`, resolve to a file inside `root`? */
function resolves(root, from, value) {
  const path = value.split('#')[0].split('?')[0];
  if (path === '') return true;
  const target = path.startsWith('/') ? resolve(root, '.' + path) : resolve(dirname(from), path);
  const file = existsSync(target) && statSync(target).isDirectory() ? join(target, 'index.html') : target;
  return (file === root || file.startsWith(root + '/')) && existsSync(file) && statSync(file).isFile();
}

export function moduleSpecifiers(code) {
  const out = [];
  const patterns = [
    /\b(?:import|export)\s[^'"`;()]*?\bfrom\s*(['"])([^'"\n]+)\1/g,
    /\bimport\s*(['"])([^'"\n]+)\1/g,
    /\bimport\s*\(\s*(['"])([^'"\n]+)\1\s*\)/g,
  ];
  for (const re of patterns) for (const m of code.matchAll(re)) out.push(m[2]);
  return out;
}

export function cssReferences(css) {
  const out = [];
  for (const m of css.matchAll(/url\(\s*(['"]?)([^'")]+)\1\s*\)/gi)) out.push(m[2].trim());
  for (const m of css.matchAll(/@import\s+(['"])([^'"]+)\1/gi)) out.push(m[2]);
  return out;
}

export function checkDist(rootArg) {
  const root = resolve(rootArg);
  const files = walk(root);
  const problems = [];
  const rel = (f) => relative(root, f);
  for (const f of files) {
    const ext = extname(f);
    if (TEXT.has(ext) && readFileSync(f, 'utf8').includes('__BUILD__')) problems.push(`${rel(f)}: __BUILD__ token left`);
  }
  const attr = /<(a|link|script|img|source)\b[^>]*?\s(href|src)\s*=\s*("([^"]*)"|'([^']*)')/gi;
  for (const f of files) {
    const ext = extname(f);
    const text = TEXT.has(ext) ? readFileSync(f, 'utf8') : '';
    let refs = [];
    if (ext === '.html') refs = [...text.matchAll(attr)].map((m) => m[4] ?? m[5]);
    else if (ext === '.js' || ext === '.mjs') refs = moduleSpecifiers(text).filter((v) => /^(\.\.?\/|\/)/.test(v));
    else if (ext === '.css') refs = cssReferences(text);
    for (const v of refs) {
      if (skip(v)) continue;
      if (!resolves(root, f, v)) problems.push(`${rel(f)}: ${v} does not resolve to a file in dist-site`);
    }
  }
  const landingPath = join(root, 'index.html');
  if (!existsSync(landingPath)) problems.push('index.html is missing');
  else {
    const landing = readFileSync(landingPath, 'utf8');
    const button = landing.match(/<a\b[^>]*\bclass="[^"]*\bbutton\b[^"]*"[^>]*>[^<]*<\/a>/i);
    if (!button) problems.push('index.html: no button link');
    else if (!/\shref="analyze\/"/.test(button[0])) problems.push(`index.html: the button does not point to analyze/: ${button[0]}`);
  }
  if (!existsSync(join(root, 'analyze', 'index.html'))) problems.push('analyze/index.html is missing');
  return problems;
}
