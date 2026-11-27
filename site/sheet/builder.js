// Table model of the sample sheet builder (site/sheet/). The pairing of file names into samples is done by the Rust rules
// (qrnamap-core fastq_names, through group_fastq_names of the web package); this file holds the edits the user makes after
// that: rename, covariate columns, columns from name parts, pair two single-end rows, split a pair, swap mates, and the
// sheet text. No imports: the unit tests load this file in Node too.

let lastKey = 0;
const nextKey = () => ++lastKey;

/** The file name of a path: the part after the last `/` or `\`. */
export const baseName = (path) => String(path).split(/[\\/]/).pop();

/** Table rows from the `samples` of group_fastq_names. */
export function rowsFromGrouping(samples) {
  return samples.map((s) => ({
    key: nextKey(),
    name: s.name,
    fromName: s.fromName,
    fastq1: [...s.fastq1],
    fastq2: [...s.fastq2],
    rule: s.rule,
    problems: [...s.problems],
    notes: [...s.notes],
    values: {},
  }));
}

/** A new covariate column. */
export const newColumn = (name) => ({ key: nextKey(), name });

/**
 * Parts of each sample name split at `_`, `-` and `.` (`WT_rep1_30C` -> WT, rep1, 30C). Returns
 * {count, parts}: `count` is the largest number of parts, and shorter names get empty parts at the end.
 */
export function nameParts(names) {
  const parts = names.map((n) => String(n).split(/[_.-]+/).filter(Boolean));
  const count = Math.max(0, ...parts.map((p) => p.length));
  return { count, parts: parts.map((p) => [...p, ...Array(count - p.length).fill('')]) };
}

/**
 * Add one covariate column per name part that differs between samples (a part that is the same in every row says nothing).
 * Returns the new columns. `rows` are changed in place.
 */
export function columnsFromNameParts(rows, columns) {
  const { count, parts } = nameParts(rows.map((r) => r.name));
  const added = [];
  for (let i = 0; i < count; i++) {
    const values = parts.map((p) => p[i]);
    if (new Set(values).size < 2) continue;
    let name = `part${i + 1}`;
    while (columns.some((c) => c.name === name) || added.some((c) => c.name === name)) name += '_';
    const col = newColumn(name);
    rows.forEach((r, j) => { r.values[col.key] = values[j]; });
    added.push(col);
  }
  return added;
}

/** Swap mate 1 and mate 2 of a paired row (a single-end row is left as it is). */
export function swapMates(row) {
  if (!row.fastq2.length) return false;
  [row.fastq1, row.fastq2] = [row.fastq2, row.fastq1];
  return true;
}

/** Split a paired row into two single-end rows (`name-1`, `name-2`). Returns the new rows, or null for a single-end row. */
export function splitRow(row) {
  if (!row.fastq1.length || !row.fastq2.length) return null;
  const make = (files, n) => ({
    ...row, key: nextKey(), name: `${row.name}-${n}`, fastq1: files, fastq2: [], rule: 'single',
    problems: [], notes: [`split from the pair "${row.name}"`], values: { ...row.values },
  });
  return [make(row.fastq1, 1), make(row.fastq2, 2)];
}

/**
 * Pair two single-end rows as mate 1 (`a`) and mate 2 (`b`). The pair takes the name and values of `a`.
 * Returns the new row, or an error message.
 */
export function pairRows(a, b) {
  if (a.fastq2.length || b.fastq2.length) return 'Pair takes two single-end rows (rows with no mate 2 file).';
  if (a.fastq1.length !== b.fastq1.length) return `"${a.name}" has ${a.fastq1.length} files and "${b.name}" has ${b.fastq1.length}: mates need the same number of lanes.`;
  return {
    ...a, key: nextKey(), fastq1: a.fastq1, fastq2: b.fastq1, rule: 'manual', problems: [],
    notes: [`paired by hand: mate 2 from "${b.name}"`], values: { ...a.values },
  };
}

/** The path written in the sheet: `prefix` + the file name, for a file given without a folder. */
export function sheetPath(path, prefix) {
  const p = String(prefix || '').trim();
  if (!p || baseName(path) !== path) return path;
  return p.endsWith('/') || p.endsWith('\\') ? p + path : `${p}/${path}`;
}

/**
 * The sheet text (the command-line format: sample, fastq1, fastq2, covariates; lanes as comma-separated lists). Every
 * column is written, also an empty one: the check then names the empty cells to fill.
 */
export function sheetText(columns, rows, prefix = '') {
  const clean = (x) => String(x ?? '').replace(/[\t\r\n]+/g, ' ').trim();
  const head = ['sample', 'fastq1', 'fastq2', ...columns.map((c) => clean(c.name))];
  const lines = rows.map((r) => [
    clean(r.name),
    r.fastq1.map((f) => sheetPath(f, prefix)).join(','),
    r.fastq2.map((f) => sheetPath(f, prefix)).join(','),
    ...columns.map((c) => clean(r.values[c.key])),
  ].join('\t'));
  return [head.join('\t'), ...lines].join('\n') + '\n';
}

/** Messages about rows that the page (/analyze/) cannot take as they are; the command-line tool can. */
export function analyzePageNotes(rows) {
  const lanes = rows.filter((r) => r.fastq1.length > 1 || r.fastq2.length > 1).map((r) => r.name);
  return lanes.length
    ? [`${lanes.length} sample(s) have several lanes (${lanes.slice(0, 3).join(', ')}${lanes.length > 3 ? ', …' : ''}). ` +
      'The command-line tool takes them as they are; the Analyze page takes one file per mate, so join the lanes first (cat L001 L002 > one file).']
    : [];
}

/** A plain copy of the model for localStorage (and back). */
export const toSaved = (state) => JSON.parse(JSON.stringify(state));
export function fromSaved(saved) {
  if (!saved || !Array.isArray(saved.rows) || !Array.isArray(saved.columns) || !Array.isArray(saved.names)) return null;
  for (const x of [...saved.rows, ...saved.columns]) lastKey = Math.max(lastKey, Number(x.key) || 0);
  return saved;
}
