// Sample sheet helpers of the samples step. The sheet format is the one of the command-line tool (crates/qrnamap-core/src/sheet.rs);
// `validate_sheet` of qrnamap-wasm checks it (the same parser). This file adds the rules of the page:
// - one file per mate (no comma-separated lanes) in this version;
// - files are matched by file name only (the folder part of a sheet path is ignored), so each file name can be used once,
//   and two picked files cannot have the same name.
// No imports: the tests load this file in Node too.

/** The file name of a path: the part after the last `/` or `\`. */
export const baseName = (path) => String(path).split(/[\\/]/).pop();

/**
 * The page rules for the rows of a valid sheet (`rows` of validate_sheet: {sample, fastq1: [], fastq2: [], covariates}).
 * Returns a list of messages (empty when the rows can be used).
 */
export function pageRules(rows) {
  const errors = [];
  const used = new Map(); // file name -> where it is used
  for (const row of rows) {
    for (const column of ['fastq1', 'fastq2']) {
      const cell = row[column];
      if (cell.length > 1) {
        errors.push(`sample "${row.sample}": the ${column} cell lists ${cell.length} files (lanes). This version of the page takes one file per mate. ` +
          'Join the lanes into one file first (gzip files can be joined with cat), or use the command-line tool.');
        continue;
      }
      if (!cell.length) continue;
      const name = baseName(cell[0]);
      const where = `sample "${row.sample}" ${column}`;
      if (used.has(name)) {
        errors.push(`the file name "${name}" is used twice (${used.get(name)} and ${where}). The page matches files by name only, so each file name can be used once.`);
      } else {
        used.set(name, where);
      }
    }
  }
  return errors;
}

/** Picked files by name. Two files with the same name are an error, because the page could not tell them apart. */
export function filesByName(files) {
  const byName = new Map();
  const errors = [];
  for (const file of files) {
    if (byName.has(file.name)) errors.push(`two picked files have the same name "${file.name}". The page matches files by name only: pick only one of them.`);
    else byName.set(file.name, file);
  }
  return { byName, errors };
}

/**
 * A sheet in the command-line format from the table of the page. `columns`: covariate column headers (a header can end
 * with `:numeric`); a column with no value in any row is left out. `rows`: {name, values: [text per column], fastq1, fastq2}
 * (file names; fastq2 empty for a single-end sample).
 */
export function sheetText(columns, rows) {
  const keep = columns.map((_, i) => rows.some((r) => String(r.values[i] ?? '').trim() !== ''));
  const head = ['sample', 'fastq1', 'fastq2', ...columns.filter((_, i) => keep[i])];
  const lines = rows.map((r) => [r.name, r.fastq1, r.fastq2, ...r.values.filter((_, i) => keep[i])].map((x) => String(x ?? '').trim()).join('\t'));
  return [head.join('\t'), ...lines].join('\n') + '\n';
}
