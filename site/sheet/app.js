// The sample sheet builder page (site/sheet/index.html). File names come from picked files, a dropped folder or pasted text;
// group_fastq_names (Rust rules, web package) pairs them into samples; the table edits of builder.js follow; validate_sheet
// (the command-line parser) checks the sheet before it is downloaded. Only file names are read, never file contents.
// The draft (names, columns, rows, folder) is kept in localStorage so a reload does not lose the work.
import init, { group_fastq_names, validate_sheet } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import {
  baseName, rowsFromGrouping, newColumn, columnsFromNameParts, swapMates, splitRow, pairRows, sheetText,
  analyzePageNotes, toSaved, fromSaved,
} from './builder.js?v=__BUILD__';

const STORE = 'qrnamap.sheetBuilder';
const $ = (id) => document.getElementById(id);
const fmt = (n) => Number(n).toLocaleString('en-US');

let state = { names: [], numericMates: true, columns: [newColumn('condition')], rows: [], skipped: [], prefix: '' };
const selected = new Set(); // row keys

function save() {
  try { localStorage.setItem(STORE, JSON.stringify(toSaved(state))); } catch { /* private window: no draft */ }
}
function load() {
  try { return fromSaved(JSON.parse(localStorage.getItem(STORE) || 'null')); } catch { return null; }
}

const el = (tag, props = {}, children = []) => {
  const e = Object.assign(document.createElement(tag), props);
  e.append(...children);
  return e;
};
const showList = (id, messages) => {
  $(id).replaceChildren(...messages.map((m) => el('li', { textContent: m })));
  $(id).hidden = messages.length === 0;
};

// ------------------------------------------------------------------ names -> rows
/** Group all names again. Column values of a sample with the same name as before are kept. */
function regroup() {
  const result = JSON.parse(group_fastq_names(JSON.stringify({ names: state.names, numericMates: state.numericMates })));
  if (!result.ok) { showList('action-error', [result.error]); return; }
  const before = new Map(state.rows.map((r) => [r.name, r.values]));
  state.rows = rowsFromGrouping(result.samples);
  for (const r of state.rows) if (before.has(r.name)) r.values = { ...before.get(r.name) };
  state.skipped = result.skipped;
  selected.clear();
  renderAll();
}

function addNames(names) {
  const known = new Set(state.names);
  const fresh = names.map((n) => String(n).trim()).filter((n) => n && !known.has(n) && (known.add(n), true));
  state.names.push(...fresh);
  $('names-status').textContent = `${fmt(fresh.length)} name(s) added, ${fmt(state.names.length)} in all.` +
    (names.length > fresh.length ? ` ${fmt(names.length - fresh.length)} already listed or empty.` : '');
  regroup();
}

/** File names of a drop: files, and the files of dropped folders (all levels), as folder-relative paths. */
async function droppedNames(items) {
  const names = [];
  const walk = async (entry, path) => {
    if (entry.isFile) { names.push(path + entry.name); return; }
    if (!entry.isDirectory) return;
    const reader = entry.createReader();
    for (;;) {
      const batch = await new Promise((res, rej) => reader.readEntries(res, rej));
      if (!batch.length) break;
      for (const e of batch) await walk(e, `${path}${entry.name}/`);
    }
  };
  const entries = [...items].map((i) => i.webkitGetAsEntry?.()).filter(Boolean);
  for (const e of entries) await walk(e, '');
  return names;
}

// ------------------------------------------------------------------ render
function renderSummary() {
  const rows = state.rows;
  const bad = rows.filter((r) => r.problems.length).length;
  const paired = rows.filter((r) => !r.problems.length && r.fastq2.length).length;
  const single = rows.length - bad - paired;
  $('group-summary').textContent = `${fmt(rows.length)} sample(s): ${fmt(paired)} paired-end, ${fmt(single)} single-end` +
    (bad ? `, ${fmt(bad)} with a problem (marked below).` : '.');
  $('skipped-box').hidden = !state.skipped.length;
  $('skipped-summary').textContent = `${fmt(state.skipped.length)} file(s) left out`;
  $('skipped-list').replaceChildren(...state.skipped.map((s) => el('li', {}, [el('code', { textContent: s.path }), `: ${s.reason}`])));
}

function fileCell(files) {
  if (!files.length) return el('td', { className: 'files small', textContent: '—' });
  return el('td', { className: 'files' }, files.map((f) => el('div', { textContent: baseName(f), title: f })));
}

function renderTable() {
  const allChecked = state.rows.length > 0 && state.rows.every((r) => selected.has(r.key));
  const head = el('tr', {}, [
    el('th', {}, [el('input', { type: 'checkbox', id: 'select-all', checked: allChecked, ariaLabel: 'Select all rows' })]),
    el('th', { textContent: 'Sample' }),
  ]);
  for (const c of state.columns) {
    const name = el('input', { type: 'text', className: 'col-name', value: c.name, size: 10, ariaLabel: `Name of column ${c.name}`, spellcheck: false });
    name.dataset.col = c.key;
    const rm = el('button', { type: 'button', className: 'remove-col', textContent: '×', title: `Remove column ${c.name}`, ariaLabel: `Remove column ${c.name}` });
    rm.dataset.col = c.key;
    head.append(el('th', {}, [name, rm]));
  }
  head.append(el('th', { textContent: 'Mate 1 (or single-end)' }), el('th', { textContent: 'Mate 2' }), el('th', { textContent: 'Notes' }));

  const body = state.rows.map((r, i) => {
    const tr = el('tr', { className: r.problems.length ? 'has-problem' : '' });
    tr.dataset.row = r.key;
    const box = el('input', { type: 'checkbox', className: 'row-select', checked: selected.has(r.key), ariaLabel: `Select row ${i + 1} (${r.name})` });
    tr.append(el('td', {}, [box]));
    tr.append(el('td', {}, [el('input', { type: 'text', className: 'sample-name', value: r.name, size: 14, ariaLabel: `Sample name, row ${i + 1}`, spellcheck: false })]));
    for (const c of state.columns) {
      const input = el('input', { type: 'text', className: 'cov', value: r.values[c.key] ?? '', size: 8, ariaLabel: `${c.name}, row ${i + 1}`, spellcheck: false });
      input.dataset.col = c.key;
      tr.append(el('td', {}, [input]));
    }
    tr.append(fileCell(r.fastq1), fileCell(r.fastq2));
    const notes = [...r.problems.map((p) => el('div', { className: 'problem', textContent: p })), ...r.notes.map((n) => el('div', { className: 'small', textContent: n }))];
    tr.append(el('td', { className: 'notes' }, notes));
    return tr;
  });
  $('builder-table').replaceChildren(el('thead', {}, [head]), el('tbody', {}, body));
  renderSelection();
}

function renderSelection() {
  $('selected-count').textContent = `${selected.size} selected`;
  const opts = state.columns.map((c) => el('option', { value: c.key, textContent: c.name || '(no name)' }));
  const keep = $('bulk-col').value;
  $('bulk-col').replaceChildren(...opts);
  if (state.columns.some((c) => String(c.key) === keep)) $('bulk-col').value = keep;
  const sa = $('select-all');
  if (sa) sa.checked = state.rows.length > 0 && state.rows.every((r) => selected.has(r.key));
}

function renderSheet() {
  const text = sheetText(state.columns, state.rows, state.prefix);
  $('sheet-preview').textContent = text;
  const messages = [];
  for (const r of state.rows) for (const p of r.problems) messages.push(`sample "${r.name}": ${p}`);
  if (!state.rows.length) messages.push('No samples yet.');
  else {
    const check = JSON.parse(validate_sheet(text));
    if (!check.ok) messages.push(...check.errors);
  }
  showList('sheet-check', messages);
  $('sheet-ok').hidden = messages.length > 0;
  $('sheet-ok').textContent = `The sheet is valid: ${fmt(state.rows.length)} samples, ${fmt(state.columns.length)} covariate column(s).`;
  showList('page-notes', analyzePageNotes(state.rows));
  save();
}

function renderAll() {
  const any = state.rows.length > 0 || state.skipped.length > 0;
  $('sheet-table-step').hidden = !any;
  $('sheet-out').hidden = !any;
  showList('action-error', []);
  renderSummary();
  renderTable();
  renderSheet();
}

// ------------------------------------------------------------------ events
function wire() {
  $('pick-files').addEventListener('change', ({ target }) => { addNames([...target.files].map((f) => f.name)); target.value = ''; });
  $('pick-folder').addEventListener('change', ({ target }) => {
    addNames([...target.files].map((f) => f.webkitRelativePath || f.name));
    target.value = '';
  });
  $('btn-use-pasted').addEventListener('click', () => {
    addNames($('paste-names').value.split(/\r?\n/));
    $('paste-names').value = '';
  });
  const zone = $('drop-zone');
  for (const ev of ['dragenter', 'dragover']) zone.addEventListener(ev, (e) => { e.preventDefault(); zone.classList.add('over'); });
  zone.addEventListener('dragleave', () => zone.classList.remove('over'));
  zone.addEventListener('drop', async (e) => {
    e.preventDefault();
    zone.classList.remove('over');
    const names = e.dataTransfer.items ? await droppedNames(e.dataTransfer.items) : [...e.dataTransfer.files].map((f) => f.name);
    addNames(names);
  });
  $('numeric-mates').addEventListener('change', ({ target }) => { state.numericMates = target.checked; regroup(); });
  $('btn-start-over').addEventListener('click', () => {
    state = { names: [], numericMates: $('numeric-mates').checked, columns: [newColumn('condition')], rows: [], skipped: [], prefix: '' };
    $('path-prefix').value = '';
    $('names-status').textContent = 'Cleared.';
    selected.clear();
    renderAll();
  });

  const table = $('builder-table');
  const rowOf = (t) => state.rows.find((r) => String(r.key) === t.closest('tr')?.dataset.row);
  table.addEventListener('input', ({ target }) => {
    if (target.classList.contains('col-name')) {
      state.columns.find((c) => String(c.key) === target.dataset.col).name = target.value.trim();
      renderSelection();
    } else {
      const r = rowOf(target);
      if (!r) return;
      if (target.classList.contains('sample-name')) r.name = target.value.trim();
      if (target.classList.contains('cov')) r.values[target.dataset.col] = target.value.trim();
    }
    renderSheet();
  });
  table.addEventListener('change', ({ target }) => {
    if (target.id === 'select-all') {
      for (const r of state.rows) target.checked ? selected.add(r.key) : selected.delete(r.key);
      table.querySelectorAll('.row-select').forEach((b) => { b.checked = target.checked; });
    } else if (target.classList.contains('row-select')) {
      const r = rowOf(target);
      target.checked ? selected.add(r.key) : selected.delete(r.key);
    }
    renderSelection();
  });
  table.addEventListener('click', ({ target }) => {
    if (!target.classList.contains('remove-col')) return;
    const key = target.dataset.col;
    state.columns = state.columns.filter((c) => String(c.key) !== key);
    for (const r of state.rows) delete r.values[key];
    renderTable();
    renderSheet();
  });

  $('btn-add-col').addEventListener('click', () => {
    const name = $('new-col-name').value.trim() || `column${state.columns.length + 1}`;
    state.columns.push(newColumn(name));
    $('new-col-name').value = '';
    renderTable();
    renderSheet();
  });
  $('btn-name-parts').addEventListener('click', () => {
    const added = columnsFromNameParts(state.rows, state.columns);
    if (!added.length) { showList('action-error', ['The sample names have no part (between _ - .) that differs between samples.']); return; }
    state.columns.push(...added);
    renderTable();
    renderSheet();
  });

  const chosen = () => state.rows.filter((r) => selected.has(r.key));
  const fail = (m) => showList('action-error', [m]);
  const done = () => { showList('action-error', []); renderSummary(); renderTable(); renderSheet(); };
  $('btn-bulk-set').addEventListener('click', () => {
    const rows = chosen();
    if (!rows.length) return fail('Select rows first (the boxes on the left).');
    const key = $('bulk-col').value;
    if (!key) return fail('Add a column first.');
    for (const r of rows) r.values[key] = $('bulk-value').value.trim();
    done();
  });
  $('btn-pair').addEventListener('click', () => {
    const rows = chosen();
    if (rows.length !== 2) return fail('Select exactly two single-end rows to pair (the upper row becomes mate 1).');
    const made = pairRows(rows[0], rows[1]);
    if (typeof made === 'string') return fail(made);
    state.rows = state.rows.filter((r) => r !== rows[1]);
    state.rows[state.rows.indexOf(rows[0])] = made;
    selected.clear();
    done();
  });
  $('btn-split').addEventListener('click', () => {
    const rows = chosen().filter((r) => r.fastq1.length && r.fastq2.length);
    if (!rows.length) return fail('Select one or more paired rows to split.');
    for (const r of rows) state.rows.splice(state.rows.indexOf(r), 1, ...splitRow(r));
    selected.clear();
    done();
  });
  $('btn-swap').addEventListener('click', () => {
    const rows = chosen();
    if (!rows.length || !rows.map(swapMates).some(Boolean)) return fail('Select one or more paired rows to swap.');
    done();
  });
  $('btn-remove').addEventListener('click', () => {
    if (!selected.size) return fail('Select rows first (the boxes on the left).');
    state.rows = state.rows.filter((r) => !selected.has(r.key));
    selected.clear();
    done();
  });

  $('path-prefix').addEventListener('input', ({ target }) => { state.prefix = target.value; renderSheet(); });
  $('btn-download').addEventListener('click', () => {
    const blob = new Blob([sheetText(state.columns, state.rows, state.prefix)], { type: 'text/tab-separated-values' });
    const a = el('a', { href: URL.createObjectURL(blob), download: 'samples.tsv' });
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  });
  $('btn-copy').addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(sheetText(state.columns, state.rows, state.prefix));
      $('copy-status').textContent = 'Copied.';
    } catch {
      $('copy-status').textContent = 'The browser did not allow copying: select the preview text instead.';
    }
  });
}

await init({ module_or_path: new URL('../pkg-web/qrnamap_wasm_bg.wasm?v=__BUILD__', import.meta.url) });
wire();
const saved = load();
if (saved) {
  state = { ...state, ...saved };
  $('numeric-mates').checked = state.numericMates;
  $('path-prefix').value = state.prefix || '';
  if (state.names.length) $('names-status').textContent = `${fmt(state.names.length)} name(s) from your last visit. "Start over" clears them.`;
}
renderAll();
document.body.dataset.ready = 'true';
