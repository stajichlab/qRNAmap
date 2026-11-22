// Link test for the built site: node check-links.mjs <absolute path of dist-site>
// Parses every HTML file in the directory. Checks that each local href/src resolves to a file,
// that the "Analyze reads" button points to analyze/, and that no __BUILD__ token is left in any text file.
// Needs no packages.
import { readFileSync, readdirSync, statSync, existsSync } from 'node:fs';
import { join, resolve, dirname, extname, relative } from 'node:path';

const root = resolve(process.argv[2] ?? '');
if (!process.argv[2] || !statSync(root, { throwIfNoEntry: false })?.isDirectory()) {
  console.error('usage: node check-links.mjs <absolute path of dist-site>');
  process.exit(2);
}
const walk = (dir) => readdirSync(dir).flatMap((n) => {
  const p = join(dir, n);
  return statSync(p).isDirectory() ? walk(p) : [p];
});
const files = walk(root);
const problems = [];
const TEXT = new Set(['.html', '.js', '.mjs', '.css', '.json', '.md', '.txt']);

for (const f of files) {
  if (TEXT.has(extname(f)) && readFileSync(f, 'utf8').includes('__BUILD__')) problems.push(`${relative(root, f)}: __BUILD__ token left`);
}

const attr = /<(a|link|script|img)\b[^>]*?\s(href|src)\s*=\s*("([^"]*)"|'([^']*)')/gi;
let checked = 0;
for (const f of files.filter((x) => x.endsWith('.html'))) {
  const html = readFileSync(f, 'utf8');
  for (const m of html.matchAll(attr)) {
    const value = m[4] ?? m[5];
    if (/^(https?:|mailto:|data:|#)/.test(value) || value === '') continue;
    const path = value.split('#')[0].split('?')[0];
    const target = resolve(path.startsWith('/') ? root : dirname(f), '.' + (path.startsWith('/') ? path : '/' + path));
    const file = existsSync(target) && statSync(target).isDirectory() ? join(target, 'index.html') : target;
    checked++;
    if (!file.startsWith(root) || !existsSync(file)) problems.push(`${relative(root, f)}: ${value} does not resolve to a file in dist-site`);
  }
}

const landing = readFileSync(join(root, 'index.html'), 'utf8');
const button = landing.match(/<a\b[^>]*\bclass="[^"]*\bbutton\b[^"]*"[^>]*>\s*Analyze reads\s*<\/a>/i);
if (!button) problems.push('index.html: no "Analyze reads" button');
else if (!/\shref="analyze\/"/.test(button[0])) problems.push(`index.html: the button does not point to analyze/: ${button[0]}`);
if (!existsSync(join(root, 'analyze', 'index.html'))) problems.push('analyze/index.html is missing');

if (problems.length) {
  console.error(problems.join('\n'));
  process.exit(1);
}
console.log(`check-links: ${checked} local links resolve, the button points to analyze/, no __BUILD__ token (${files.length} files)`);
