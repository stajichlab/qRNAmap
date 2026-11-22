// Link test for the built site: node check-links.mjs <absolute path of dist-site>
// The checks are in links-lib.mjs. Exit 1 when a problem is found.
import { statSync } from 'node:fs';
import { checkDist } from './links-lib.mjs';

const root = process.argv[2];
if (!root || !statSync(root, { throwIfNoEntry: false })?.isDirectory()) {
  console.error('usage: node check-links.mjs <absolute path of dist-site>');
  process.exit(2);
}
const problems = checkDist(root);
if (problems.length) { console.error(problems.join('\n')); process.exit(1); }
console.log('check-links: all local links, module imports and CSS urls resolve; the button points to analyze/; no __BUILD__ token');
