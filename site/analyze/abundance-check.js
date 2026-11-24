// Check one kallisto abundance.tsv before the page uses it. Reason: in the Task 3 measurements Firefox wrote a damaged
// abundance.tsv in 4 of 130 kallisto runs, with no error (docs/browser/page-v1-measurements.md, "Run-to-run output of
// kallisto"). The damage seen: a number token written more than once (`9` became `999`; `204.998` became
// `204.998204.998204.998`). The check there would have caught both cases:
//   1. the header is the kallisto header, and every data line has a name and 4 numbers;
//   2. the number of data lines equals the number of transcripts in the index;
//   3. the sum of est_counts equals n_pseudoaligned of run_info.json, within the rounding of the printed numbers.
// Rule 3 tolerance: kallisto prints each est_counts with 6 significant digits, so each printed value can be off by half of
// its last digit. The tolerance is the sum of these halves, plus 1. Not caught: a repeated token that changes only digits
// far after the decimal point (it changes the sum by less than the tolerance). No imports: this file runs in the page,
// in the worker and in Node (unit tests).

const HEADER = 'target_id\tlength\teff_length\test_counts\ttpm';
const NUMBER = /^[+-]?(\d+)(?:\.(\d+))?(?:[eE]([+-]?\d+))?$/;
// tpm is 0/0 when no read pseudoaligns; kallisto then prints nan. Only that column may hold it.
const NAN = /^[+-]?nan$/i;

/** Half of the last printed digit of a number text (for example 0.0005 for "1.234"). */
function halfDigit(text) {
  const m = NUMBER.exec(text);
  const decimals = m[2] ? m[2].length : 0;
  const exponent = m[3] ? Number(m[3]) : 0;
  return 0.5 * 10 ** (exponent - decimals);
}

/**
 * Check the text of an abundance.tsv. Returns a list of problems in plain words (empty when the file passes).
 * `nTranscripts`: the transcripts in the index. `nPseudoaligned`: n_pseudoaligned from run_info.json.
 */
export function checkAbundance(text, { nTranscripts, nPseudoaligned }) {
  const problems = [];
  const lines = text.split('\n');
  if (lines[lines.length - 1] === '') lines.pop();
  if (lines[0] !== HEADER) problems.push('the first line is not the kallisto header');
  const rows = lines.slice(1);
  if (rows.length !== nTranscripts) problems.push(`it has ${rows.length} transcript lines, the index has ${nTranscripts} transcripts`);
  let sum = 0;
  let tolerance = 1;
  let bad = 0;
  for (const [i, row] of rows.entries()) {
    const f = row.split('\t');
    const numbers = f.slice(1);
    const ok = f.length === 5 && f[0] !== '' && numbers.every((x, k) => NUMBER.test(x) || (k === 3 && NAN.test(x)));
    if (!ok) {
      if (bad++ === 0) problems.push(`line ${i + 2} is not a name and 4 numbers: "${row.slice(0, 80)}"`);
      continue;
    }
    sum += Number(f[3]);
    tolerance += halfDigit(f[3]);
  }
  if (bad > 1) problems.push(`${bad} lines in all are not a name and 4 numbers`);
  if (!Number.isFinite(nPseudoaligned)) problems.push('run_info.json has no n_pseudoaligned number');
  else if (bad === 0 && Math.abs(sum - nPseudoaligned) > tolerance) {
    problems.push(`the est_counts add up to ${sum.toFixed(3)}, but run_info.json says ${nPseudoaligned} reads pseudoaligned`);
  }
  return problems;
}
