// Step 2 of /analyze/: the sample table (sample name, covariate columns, one mate 1 and one mate 2 file per row), the
// upload of a command-line sample sheet with file-name matching, the checks before a run, and the run (quant.js).
// The table is the one place the samples live: a loaded sheet fills it, and a run reads it.
// Checks: the sheet that the table gives is checked with validate_sheet (the command-line parser and reserved names)
// and with the page rules of sheet.js. The page does not refuse a file by its size: it reads only the first N pairs.
import { validate_sheet } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { baseName, pageRules, filesByName, sheetText } from './sheet.js?v=__BUILD__';
import { createQuant, MAX_PAIRS, previewLabel } from './quant.js?v=__BUILD__';
import { isTauri, checkNative, pickFile, pickFiles, readLocalSheet, FASTQ_EXTS } from './tauri-bridge.js?v=__BUILD__';

const fmt = (n) => Number(n).toLocaleString('en-US');
// Measurement option (page-quant-memory.yml): `?preview_pairs=N` replaces the read-length limit of preview-limit.js with N for every
// sample, so a candidate N can be measured through the shipped page path. It is capped at the hard limit (200,000) by the worker.
const askedOverride = Math.floor(Number(new URLSearchParams(location.search).get('preview_pairs')));
const LIMIT_OVERRIDE = Number.isFinite(askedOverride) && askedOverride >= 1 ? Math.min(askedOverride, 200000) : 0;
const ROUTE = 'Use the command-line tool instead (README, "Install and build"): qrnamap quant --samples samples.tsv ...';

/**
 * Set up step 2. `$` finds an element by id. `app`: window.qrnamapApp ({state, getWorker, nextId}).
 * `shell`: { markDone(name), resetFrom(name) } of app.js.
 */
export function createSamplesStep({ $, app, shell }) {
  let nextKey = 0;
  const key = () => ++nextKey;
  const model = { columns: [{ key: key(), name: 'condition' }], rows: [] };
  let picked = new Map(); // files picked for a sheet, by name
  let runRows = []; // the table rows of the current run, in run order

  const newRow = (fields = {}) => ({ key: key(), name: '', values: {}, files: { 1: null, 2: null }, expect: { 1: null, 2: null }, ...fields });
  model.rows.push(newRow());

  const showList = (id, messages) => {
    const box = $(id);
    box.replaceChildren(...messages.map((m) => Object.assign(document.createElement('li'), { textContent: m })));
    box.hidden = messages.length === 0;
  };
  const sheetError = (messages) => showList('sheet-error', messages);

  // ------------------------------------------------------------------ the table
  const el = (tag, { dataset = {}, ...props } = {}, children = []) => {
    const e = Object.assign(document.createElement(tag), props);
    for (const [k, v] of Object.entries(dataset)) e.dataset[k] = v;
    e.append(...children);
    return e;
  };

  function mateCell(row, mate) {
    const file = row.files[mate];
    const expect = row.expect[mate];
    let text = '';
    let missing = false;
    if (file) {
      if (row.exists && row.exists[mate] && row.exists[mate].includes(false)) {
        text = `missing on disk: ${file.path || file.name}`;
        missing = true;
      } else if (file.path) {
        text = `${file.name} (local disk)`;
      } else {
        text = expect && file.name === expect ? `${file.name} (matched to the sheet)` : file.name;
      }
    }
    else if (row.loading && row.loading[mate]) text = row.loading[mate];
    else if (expect) { text = `missing: the sheet names "${expect}". Pick that file.`; missing = true; }
    else if (mate === 2) text = 'none (single-end)';
    const input = el('input', {
      type: 'file',
      className: 'mate-file',
      ariaLabel: `Mate ${mate} file of sample row ${model.rows.indexOf(row) + 1}`,
      accept: '.fastq,.fq,.gz,.fastq.gz,application/gzip,application/x-gzip,text/plain,*/*',
    });
    input.dataset.mate = mate;
    // The input is visually hidden inside a label styled as a button (native file inputs do not take styles and wrap badly
    // in a table cell). In the desktop app a click opens the native dialog instead (handled on the table), so the row gets a disk path.
    const pick = el('label', { className: 'file-button' }, [input, file ? 'Change…' : 'Choose file…']);
    const span = el('span', { className: missing ? 'mate-name missing' : 'mate-name', textContent: text });
    if (file?.path) span.title = file.path;
    return el('td', { className: 'mate', dataset: { mate: String(mate), missing: String(missing) } }, [pick, span]);
  }

  function render() {
    const head = el('tr', {}, [el('th', { textContent: 'Sample name' })]);
    for (const col of model.columns) {
      const name = el('input', { type: 'text', className: 'col-name', value: col.name, ariaLabel: 'Column name', size: 10 });
      name.dataset.col = col.key;
      const remove = el('button', { type: 'button', className: 'remove-col', textContent: 'Remove', ariaLabel: `Remove column ${col.name}` });
      remove.dataset.col = col.key;
      head.append(el('th', {}, [name, ' ', remove]));
    }
    head.append(el('th', { textContent: 'Mate 1 file' }), el('th', { textContent: 'Mate 2 file' }), el('th', { textContent: 'Status' }), el('th', { textContent: '' }));
    const body = model.rows.map((row, i) => {
      const tr = el('tr', { dataset: { row: String(row.key), missing: String(Boolean((row.expect[1] && !row.files[1]) || (row.expect[2] && !row.files[2]))) } });
      tr.append(el('td', {}, [el('input', { type: 'text', className: 'sample-name', value: row.name, ariaLabel: `Sample name, row ${i + 1}`, size: 12 })]));
      for (const col of model.columns) {
        const input = el('input', { type: 'text', className: 'cov', value: row.values[col.key] ?? '', ariaLabel: `${col.name}, row ${i + 1}`, size: 8 });
        input.dataset.col = col.key;
        tr.append(el('td', {}, [input]));
      }
      tr.append(mateCell(row, 1), mateCell(row, 2), el('td', { className: 'status', textContent: row.statusText || '' }));
      tr.append(el('td', {}, [el('button', { type: 'button', className: 'remove-row', textContent: 'Remove', ariaLabel: `Remove sample row ${i + 1}` })]));
      return tr;
    });
    $('sample-table').replaceChildren(el('thead', {}, [head]), el('tbody', {}, body));
  }

  const rowOf = (target) => model.rows.find((r) => String(r.key) === target.closest('tr')?.dataset.row);
  $('sample-table').addEventListener('input', ({ target }) => {
    if (target.classList.contains('col-name')) { model.columns.find((c) => String(c.key) === target.dataset.col).name = target.value; return; }
    const row = rowOf(target);
    if (!row) return;
    if (target.classList.contains('sample-name')) row.name = target.value;
    if (target.classList.contains('cov')) row.values[target.dataset.col] = target.value;
  });
  $('sample-table').addEventListener('change', ({ target }) => {
    if (!target.classList.contains('mate-file')) return;
    const row = rowOf(target);
    row.files[target.dataset.mate] = target.files[0] || null;
    render();
  });
  // Desktop app: a file picked through the page's file input is a browser File with no disk path, which only the
  // in-page preview (75,000 to 200,000 pairs) can read. The native dialog gives the path, so native kallisto reads every read.
  if (isTauri()) {
    $('sample-table').addEventListener('click', async (event) => {
      const input = event.target.closest('label.file-button')?.querySelector('input.mate-file')
        || (event.target.classList?.contains('mate-file') ? event.target : null);
      if (!input) return;
      event.preventDefault();
      const row = rowOf(input);
      const mate = Number(input.dataset.mate);
      const path = await pickFile(`Mate ${mate} FASTQ for ${row.name || 'this sample'}`, 'FASTQ (.fastq, .fq, .gz)', FASTQ_EXTS);
      if (!path) return;
      row.files[mate] = { name: baseName(path), path };
      row.paths = { ...(row.paths || {}), [mate]: [path] };
      if (row.exists) row.exists = { ...row.exists, [mate]: [true] };
      render();
    });
    $('sheet-file').addEventListener('click', async (event) => {
      event.preventDefault();
      const path = await pickFile('Sample sheet', 'Sample sheet (.tsv, .txt)', ['tsv', 'txt']);
      if (!path) return;
      if ($('native-sheet-path')) $('native-sheet-path').value = path;
      await loadNativeSheet(path); // resolves the FASTQ paths against the sheet's folder and checks them on disk
    });
    $('sheet-fastqs').addEventListener('click', async (event) => {
      event.preventDefault();
      const paths = await pickFiles('The FASTQ files the sheet names', 'FASTQ (.fastq, .fq, .gz)', FASTQ_EXTS);
      if (!paths.length) return;
      picked = new Map(paths.map((p) => [baseName(p), { name: baseName(p), path: p }]));
      for (const row of model.rows) {
        for (const mate of [1, 2]) {
          const f = row.expect[mate] && picked.get(row.expect[mate]);
          if (f) { row.paths = { ...(row.paths || {}), [mate]: [f.path] }; if (row.exists) row.exists = { ...row.exists, [mate]: [true] }; }
        }
      }
      sheetError([]);
      applyPicked();
      render();
    });
  }
  $('sample-table').addEventListener('click', ({ target }) => {
    if (target.classList.contains('remove-row')) {
      model.rows = model.rows.filter((r) => r !== rowOf(target));
      render();
    } else if (target.classList.contains('remove-col')) {
      model.columns = model.columns.filter((c) => String(c.key) !== target.dataset.col);
      render();
    }
  });
  $('add-sample').addEventListener('click', () => { model.rows.push(newRow()); render(); });
  $('add-column').addEventListener('click', () => { model.columns.push({ key: key(), name: '' }); render(); });

  // ------------------------------------------------------------------ sheet upload and file matching
  function applyPicked() {
    for (const row of model.rows) {
      for (const mate of [1, 2]) {
        if (row.expect[mate] && picked.has(row.expect[mate])) row.files[mate] = picked.get(row.expect[mate]);
      }
    }
  }

  $('sheet-file').addEventListener('change', async ({ target }) => {
    const file = target.files[0];
    if (!file) return;
    let text;
    try {
      text = await file.text();
    } catch (error) {
      sheetError([`The sample sheet ${file.name} could not be read (${error?.message || error}). Pick it again.`]);
      return;
    } finally {
      target.value = '';
    }
    const result = JSON.parse(validate_sheet(text));
    if (!result.ok) { sheetError([`The sample sheet ${file.name} was not loaded:`, ...result.errors]); return; }
    const rules = pageRules(result.rows);
    if (rules.length) { sheetError([`The sample sheet ${file.name} was not loaded:`, ...rules]); return; }
    sheetError([]);
    model.columns = result.columns.map(([name, type]) => ({ key: key(), name: type === 'numeric' ? `${name}:numeric` : name }));
    model.rows = result.rows.map((r) => newRow({
      name: r.sample,
      values: Object.fromEntries(model.columns.map((c, i) => [c.key, r.covariates[result.columns[i][0]] ?? ''])),
      expect: { 1: baseName(r.fastq1[0]), 2: r.fastq2.length ? baseName(r.fastq2[0]) : null },
    }));
    applyPicked();
    render();
    $('samples-panel').dataset.sheet = file.name;
  });

  $('sheet-fastqs').addEventListener('change', ({ target }) => {
    const { byName, errors } = filesByName([...target.files]);
    if (errors.length) { sheetError(['The picked FASTQ files were not used:', ...errors]); target.value = ''; return; }
    sheetError([]);
    picked = byName;
    applyPicked();
    render();
  });

  // ------------------------------------------------------------------ desktop native local sheet loader
  async function loadNativeSheet(sheetPath) {
    const statusEl = $('native-sheet-status');
    if (statusEl) statusEl.textContent = 'Reading local sample sheet...';
    sheetError([]);
    try {
      const res = await readLocalSheet(sheetPath);
      if (!res || !res.samples || !res.samples.length) {
        if (statusEl) statusEl.textContent = 'Sample sheet is empty or contains no samples.';
        return;
      }
      model.columns = res.columns.map(([name, type]) => ({
        key: key(),
        name: type === 'numeric' ? `${name}:numeric` : name,
      }));
      model.rows = res.samples.map((s) => newRow({
        name: s.sample || s.name || '',
        values: Object.fromEntries(model.columns.map((c, i) => [c.key, s.covariates[res.columns[i][0]] ?? ''])),
        expect: {
          1: s.fastq1[0] ? baseName(s.fastq1[0]) : null,
          2: s.fastq2[0] ? baseName(s.fastq2[0]) : null,
        },
        files: {
          1: s.fastq1[0] ? { name: baseName(s.fastq1[0]), path: s.fastq1[0] } : null,
          2: s.fastq2[0] ? { name: baseName(s.fastq2[0]), path: s.fastq2[0] } : null,
        },
        paths: {
          1: s.fastq1,
          2: s.fastq2,
        },
        exists: {
          1: s.fastq1_exists,
          2: s.fastq2_exists,
        },
      }));
      render();
      $('samples-panel').dataset.sheet = baseName(sheetPath);
      const missingCount = res.warnings.length;
      if (missingCount > 0) {
        if (statusEl) statusEl.textContent = `Loaded ${res.samples.length} samples (${missingCount} file warnings: see table or check paths).`;
        sheetError(res.warnings);
      } else {
        if (statusEl) statusEl.textContent = `Loaded ${res.samples.length} samples from ${baseName(sheetPath)}. All FASTQ files verified on disk!`;
      }
    } catch (err) {
      if (statusEl) statusEl.textContent = `Error: ${err?.message || err}`;
      sheetError([`Failed to read local sheet ${sheetPath}: ${err?.message || err}`]);
    }
  }

  $('btn-browse-sheet')?.addEventListener('click', async () => {
    try {
      const pickedPath = await pickFile('Select Sample Sheet', 'Sample sheet (*.txt, *.tsv)', ['txt', 'tsv']);
      if (pickedPath) {
        $('native-sheet-path').value = pickedPath;
        await loadNativeSheet(pickedPath);
      }
    } catch (err) {
      sheetError([`Error picking sheet: ${err.message || err}`]);
    }
  });

  $('btn-load-native-sheet')?.addEventListener('click', async () => {
    const p = $('native-sheet-path')?.value.trim();
    if (!p) {
      sheetError(['Please enter or browse for a sample sheet path first.']);
      return;
    }
    await loadNativeSheet(p);
  });

  let nativeReady = false; // native kallisto found, with an index-compatible version
  /** True when this run goes to native kallisto (desktop app, native kallisto usable, "Native" chosen). */
  const nativeSelected = () => isTauri() && nativeReady && Boolean(document.getElementById('mode-native')?.checked);

  function updateExecModeUI() {
    if (!isTauri()) return;
    const isNative = nativeSelected();
    const previewP = $('preview-pairs-p');
    const notice = $('native-pairs-notice');
    if (previewP) previewP.hidden = isNative;
    if (notice) notice.hidden = !isNative;
  }

  if (isTauri()) {
    checkNative().then((caps) => {
      nativeReady = Boolean(caps?.available && caps?.index_compatible);
      const banner = $('desktop-native-banner');
      if (banner) { banner.hidden = false; banner.dataset.native = String(nativeReady); }
      const opts = $('native-mode-options');
      if (opts) opts.hidden = !nativeReady;
      if (!nativeReady) { const b = $('mode-browser'); if (b) b.checked = true; }
      const info = $('desktop-native-info');
      if (info) {
        info.textContent = nativeReady
          ? `${String(caps.kallisto_version || '').replace(/^kallisto,?\s*(version\s*)?/i, '') || 'found'} at ${caps.kallisto_path} (${caps.num_cpus} CPUs)`
          : `not available: ${caps?.warning || caps?.error || 'kallisto was not found'} Counting uses the in-page preview until then.`;
      }
      updateExecModeUI();
    });
    document.querySelectorAll('input[name="quant-exec-mode"]').forEach((el) => {
      el.addEventListener('change', updateExecModeUI);
    });
  }

  // ------------------------------------------------------------------ checks before a run
  function check() {
    const errors = [];
    const rowErrors = new Map();
    const fail = (row, text) => { errors.push(text); if (row) rowErrors.set(row, [...(rowErrors.get(row) || []), text]); };
    if (!model.rows.length) fail(null, 'Add at least one sample.');
    const isNativeMode = nativeSelected();

    model.rows.forEach((row, i) => {
      const who = `Row ${i + 1}${row.name ? ` (${row.name})` : ''}`;
      for (const mate of [1, 2]) {
        if (row.files[mate]) {
          if (row.exists && row.exists[mate] && row.exists[mate].includes(false)) {
            fail(row, `${who}: mate ${mate} file does not exist on disk: ${row.paths?.[mate]?.join(', ')}.`);
          }
          if (isNativeMode && !row.files[mate].path && !row.paths?.[mate]?.length) {
            fail(row, `${who}: the mate ${mate} file (${row.files[mate].name}) is held in page memory (for example a test-set download), not on disk, so native kallisto cannot read it. ` +
              'Pick it again with "Choose file…" (it then opens the desktop file dialog), or choose "Browser preview" to count a preview of it.');
          }
          continue;
        }
        if (row.expect[mate]) fail(row, `${who}: the sheet names "${row.expect[mate]}" for mate ${mate}, but no picked file has that name. Pick the file.`);
        else if (mate === 1) fail(row, `${who}: pick the mate 1 file.`);
      }
      for (const mate of [1, 2]) {
        const f = row.files[mate];
        if (f && /[,\t\r\n]/.test(f.name)) {
          fail(row, `${who}: the mate ${mate} file name "${f.name}" has a comma (or a tab or line break). The sample sheet uses commas to separate lanes, so rename the file first.`);
        }
      }
    });
    model.columns.forEach((col) => {
      if (!col.name.trim() && model.rows.some((r) => String(r.values[col.key] ?? '').trim())) fail(null, 'A covariate column has values but no name. Give it a name or remove it.');
    });

    const maxAllowed = MAX_PAIRS; // the in-page preview (browser mode, also in the desktop app) never reads more than this
    const pairs = Number($('preview-pairs').value);
    if (!isNativeMode) {
      if (!Number.isInteger(pairs) || pairs < 1 || pairs > maxAllowed) {
        fail(null, `Read pairs per sample must be a whole number from 1 to ${fmt(maxAllowed)}. The page can use fewer: the read length sets the largest N of each sample.`);
      }
      const hasOnlyPaths = model.rows.some((r) => r.files[1] && !(r.files[1] instanceof File) && r.paths?.[1]);
      if (hasOnlyPaths) {
        fail(null, "Local disk files require Desktop Native Execution. Please select 'Native Execution (Full depth)'.");
      }
    }
    const single = model.rows.some((r) => r.files[1] && !r.files[2] && !r.expect[2]);
    const mean = Number($('frag-mean').value);
    const sd = Number($('frag-sd').value);
    if (single && !(mean > 0 && sd > 0 && $('frag-mean').value !== '' && $('frag-sd').value !== '')) {
      fail(null, 'Single-end samples (no mate 2 file) need the fragment length mean and SD (both above 0). kallisto needs them, as the command-line tool does (--fragment-length MEAN,SD).');
    }
    if (errors.length) return { errors, rowErrors };

    const text = sheetText(model.columns.map((c) => c.name.trim()), model.rows.map((r) => ({
      name: r.name, values: model.columns.map((c) => r.values[c.key] ?? ''), fastq1: r.files[1].name, fastq2: r.files[2]?.name ?? '',
    })));
    const result = JSON.parse(validate_sheet(text));
    const sheetErrors = result.ok ? pageRules(result.rows) : result.errors;
    for (const message of sheetErrors) {
      const rows = model.rows.filter((r) => r.name && (message.includes(`"${r.name}"`) || message.includes(`'${r.name}'`)));
      errors.push(message);
      for (const row of rows) rowErrors.set(row, [...(rowErrors.get(row) || []), message]);
    }
    if (errors.length) return { errors, rowErrors };
    return {
      errors, rowErrors, pairs, text,
      jobs: model.rows.map((r, i) => ({
        name: r.name,
        covariates: result.rows[i].covariates,
        mate1: r.files[1],
        mate2: r.files[2],
        paths1: r.paths?.[1] || (r.files[1]?.path ? [r.files[1].path] : null),
        paths2: r.paths?.[2] || (r.files[2]?.path ? [r.files[2].path] : []),
        fragment: r.files[2] ? null : [mean, sd],
      })),
    };
  }

  // ------------------------------------------------------------------ the run
  const statusCell = (row) => $('sample-table').querySelector(`tr[data-row="${row.key}"] .status`);

  function sampleText(s) {
    const unit = s.args?.includes('--single') ? 'reads' : 'pairs';
    switch (s.status) {
      case 'waiting': return 'Waiting.';
      case 'running': return 'Running.';
      case 'done': return `${fmt(s.processed)} ${unit} processed, ${fmt(s.pseudoaligned)} pseudoaligned (${(100 * s.fraction).toFixed(1)}%). ` +
        (s.allPairs
          ? `All ${unit} (${fmt(s.nPairs)}). Full depth.`
          : `${previewLabel(s.previewN).replace(/^p/, 'P')}. Reads of ${fmt(s.readLength)} bases: N = ${fmt(s.previewN)}.` +
            (isTauri() ? ' (In-page preview. Choose "Native" with files on disk for all reads.)' : '')) +
        `${s.reruns ? ' kallisto ran twice: the first output failed the check.' : ''}` +
        `${s.warning ? ` Warning: ${s.warning}` : ''}`;
      case 'error': return `Error: ${s.error}${['memory', 'kallisto'].includes(s.kind) ? ` ${ROUTE}` : ''}`;
      case 'cancelled': return 'Cancelled.';
      default: return 'Not run.';
    }
  }

  function showSample(sample) {
    const i = app.state.samples.indexOf(sample);
    const row = runRows[i];
    if (!row) return;
    const tr = $('sample-table').querySelector(`tr[data-row="${row.key}"]`);
    if (!tr) return;
    tr.dataset.status = sample.status;
    statusCell(row).textContent = sampleText(sample);
    showPreviewLabel();
  }

  /** The preview label is on the page whenever a counted sample exists (and only then). */
  function showPreviewLabel() {
    const results = (app.state.samples || []).some((s) => s.status === 'done');
    let label = document.getElementById('preview-label');
    if (!results) { label?.remove(); return; }
    if (!label) {
      label = el('p', { id: 'preview-label', className: 'preview-label' });
      $('samples-results').append(label);
    }
    label.textContent = `Counts: ${app.state.preview.label} of each sample. A sample with fewer pairs uses all its pairs ("all pairs").`;
  }

  const busy = (flag) => {
    $('samples-edit').disabled = flag;
    $('samples-busy').hidden = !flag;
    $('step-samples').dataset.phase = flag ? 'running' : (app.state.samplesPhase || 'edit');
  };

  const quant = createQuant(app, {
    onChange: showSample,
    onProgress: (p) => {
      $('samples-status').textContent = p.text;
      const i = app.state.samples.findIndex((s) => s.status === 'running');
      if (i >= 0 && runRows[i]) statusCell(runRows[i]).textContent = `Running: ${p.text}`;
    },
  });

  $('run-samples').addEventListener('click', async () => {
    showList('samples-errors', []);
    for (const tr of $('sample-table').querySelectorAll('tbody tr')) { delete tr.dataset.invalid; delete tr.dataset.status; tr.querySelector('.status').textContent = ''; }
    const checked = check();
    if (checked.errors.length) {
      showList('samples-errors', checked.errors);
      for (const [row, messages] of checked.rowErrors) {
        const tr = $('sample-table').querySelector(`tr[data-row="${row.key}"]`);
        tr.dataset.invalid = 'true';
        tr.querySelector('.status').textContent = messages.join(' ');
      }
      return;
    }
    if (!app.state.indexReady) { showList('samples-errors', ['The genome index is not ready. Build it in step 1 first.']); return; }
    shell.resetFrom('samples');
    runRows = [...model.rows];
    busy(true);
    try {
      await quant.run(checked.jobs, checked.pairs, LIMIT_OVERRIDE);
    } finally {
      busy(false);
    }
    if (app.state.samplesPhase === 'done' && app.state.samples.some((s) => s.status === 'done')) shell.markDone('samples');
    const counted = app.state.samples.filter((s) => s.status === 'done').length;
    $('samples-status').textContent = '';
    $('samples-summary').textContent = app.state.samplesPhase === 'done'
      ? `${counted} of ${app.state.samples.length} samples counted.`
      : app.state.samplesPhase === 'cancelled' ? 'Cancelled. Samples counted before the cancel keep their counts. Press "Count reads" to run all samples again.' : '';
  });
  $('samples-cancel').addEventListener('click', () => quant.cancel());

  // The genome changed, or the page worker was lost: the counts belong to the old index, so they are cleared.
  window.addEventListener('qrnamap:worker-replaced', () => {
    const had = (app.state.samples || []).length > 0;
    setTimeout(() => { // after quant.js has closed its request
      app.state.samples = [];
      app.state.samplesPhase = 'edit';
      runRows = [];
      for (const tr of $('sample-table').querySelectorAll('tbody tr')) { delete tr.dataset.status; tr.querySelector('.status').textContent = ''; }
      showPreviewLabel();
      $('samples-summary').textContent = '';
      $('step-samples').dataset.phase = 'edit';
      if (had) showList('samples-errors', ['The genome index was replaced or lost, so the counts were cleared. Build the index in step 1, then count the reads again.']);
    }, 0);
  });

  function setSamples({ columns = ['condition'], rows = [] }) {
    model.columns = columns.map((name, i) => ({ key: i + 1, name }));
    model.rows = rows.map((r, i) => ({
      key: i + 1,
      name: r.name,
      values: { 1: r.values?.[0] ?? r.condition ?? '' },
      files: { 1: r.files?.[1] ?? r.file1 ?? null, 2: r.files?.[2] ?? r.file2 ?? null },
      expect: { 1: r.expect?.[1] ?? '', 2: r.expect?.[2] ?? '' },
      loading: r.loading || null,
      statusText: r.statusText || '',
    }));
    render();
  }

  function updateRow(index, updates) {
    if (model.rows[index]) {
      Object.assign(model.rows[index], updates);
      render();
    }
  }

  $('preview-pairs').max = String(MAX_PAIRS);
  $('preview-pairs').value = String(MAX_PAIRS);
  app.state.samples = [];
  app.state.preview = { n: null, label: '' };
  app.state.samplesPhase = 'edit';
  render();
  $('step-samples').dataset.phase = 'edit';
  return { quant, model, setSamples, updateRow };
}
