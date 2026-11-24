// Test sets: load yeast AmpB (GSE80357) and C. neoformans H99 Lkh1 (PRJNA976751).
import { VerifyingGunzip } from '../gz.js?v=__BUILD__';
import { firstNameStem } from './fastq.js?v=__BUILD__';

/** Stream a gzipped or plain FASTQ URL, count newlines up to maxRecords * 4, and return the cut plain bytes. */
export async function streamCutFastq(url, maxRecords, signal) {
  const targetNewlines = maxRecords * 4;
  let res;
  try {
    res = await fetch(url, { signal });
  } catch (e) {
    throw new Error(`Failed to fetch ${url} from ENA: ${e.message || e}`);
  }
  if (!res.ok) {
    throw new Error(`ENA request failed for ${url}: HTTP ${res.status} ${res.statusText}`);
  }

  const reader = res.body.getReader();
  const chunks = [];
  let newlines = 0;
  let done = false;
  let isGzip = null;
  let gunzip = null;

  try {
    while (!done) {
      const { value, done: streamDone } = await reader.read();
      if (streamDone) break;
      if (!value || !value.length) continue;

      if (isGzip === null) {
        isGzip = value.length >= 2 && value[0] === 0x1f && value[1] === 0x8b;
        if (isGzip) {
          gunzip = new VerifyingGunzip((chunk) => {
            if (done) return;
            for (let i = 0; i < chunk.length; i++) {
              if (chunk[i] === 10) {
                newlines++;
                if (newlines === targetNewlines) {
                  chunks.push(chunk.slice(0, i + 1));
                  done = true;
                  return;
                }
              }
            }
            chunks.push(chunk);
          });
        }
      }

      if (isGzip) {
        gunzip.push(value, false);
      } else {
        // Plain FASTQ
        for (let i = 0; i < value.length; i++) {
          if (value[i] === 10) {
            newlines++;
            if (newlines === targetNewlines) {
              chunks.push(value.slice(0, i + 1));
              done = true;
              break;
            }
          }
        }
        if (!done) chunks.push(value);
      }
    }
  } finally {
    try { await reader.cancel(); } catch (_) {}
  }

  const totalLen = chunks.reduce((acc, c) => acc + c.length, 0);
  const out = new Uint8Array(totalLen);
  let off = 0;
  for (const c of chunks) {
    out.set(c, off);
    off += c.length;
  }
  return { bytes: out, records: Math.floor(newlines / 4) };
}

export function createTestsets({ $, app, shell }) {
  let testsets = [];

  async function init() {
    try {
      const res = await fetch('../testsets.json?v=__BUILD__');
      if (res.ok) testsets = await res.json();
    } catch (_) {
      // Fallback if testsets.json cannot be fetched directly
    }
  }

  function setStatus(text) {
    const el = $('testset-status');
    if (el) el.textContent = text;
  }

  function showErrors(errs) {
    const ul = $('testset-errors');
    if (!ul) return;
    if (!errs || !errs.length) {
      ul.hidden = true;
      ul.innerHTML = '';
      return;
    }
    ul.hidden = false;
    ul.innerHTML = errs.map((e) => `<li>${e}</li>`).join('');
  }

  async function loadTestSet(id) {
    showErrors([]);
    const ts = testsets.find((x) => x.id === id);
    if (!ts) {
      showErrors([`Unknown test set ${id}`]);
      return;
    }

    try {
      setStatus(`Loading test set "${ts.name}"...`);

      // 1. Fill preview pairs N
      const pairsInput = $('preview-pairs');
      if (pairsInput && ts.n_pairs) pairsInput.value = String(ts.n_pairs);

      // 2. Select catalog genome and initiate download/build
      setStatus(`Downloading and indexing catalog genome ${ts.catalog_accession} (${ts.organism})...`);
      if (app.loadCatalogAccession) {
        await app.loadCatalogAccession(ts.catalog_accession);
      }

      // 3. Populate sample table
      setStatus(`Setting up sample table with ${ts.samples.length} samples...`);
      const covName = ts.design?.replace(/^~\s*/, '') || 'condition';
      const rows = [];

      for (let i = 0; i < ts.samples.length; i++) {
        const s = ts.samples[i];
        const covVal = s[covName] || s.condition || s.genotype || '';
        setStatus(`Fetching reads from ENA for sample ${i + 1} of ${ts.samples.length} (${s.sample}, ${s.run})...`);

        // Check for test mock override function
        const url1 = window.__testSetUrl ? window.__testSetUrl(s.run, 1) : s.fastq1_url;
        const url2 = window.__testSetUrl ? window.__testSetUrl(s.run, 2) : s.fastq2_url;

        const [m1, m2] = await Promise.all([
          streamCutFastq(url1, ts.n_pairs),
          streamCutFastq(url2, ts.n_pairs),
        ]);

        if (m1.records !== m2.records) {
          throw new Error(`Record count mismatch between mates for ${s.sample}: mate 1 has ${m1.records}, mate 2 has ${m2.records}`);
        }

        const stem1 = firstNameStem(m1.bytes);
        const stem2 = firstNameStem(m2.bytes);
        if (stem1 && stem2 && stem1 !== stem2) {
          console.warn(`Read name mismatch on sample ${s.sample}: ${stem1} vs ${stem2}`);
        }

        const file1 = new File([m1.bytes], `${s.run}_1.fastq`, { type: 'text/plain' });
        const file2 = new File([m2.bytes], `${s.run}_2.fastq`, { type: 'text/plain' });

        rows.push({
          name: s.sample,
          values: [covVal],
          files: { 1: file1, 2: file2 },
          expect: { 1: file1.name, 2: file2.name },
        });
      }

      if (app.samplesStep?.setSamples) {
        app.samplesStep.setSamples({ columns: [covName], rows });
      }

      setStatus(`All samples fetched from ENA. Counting reads...`);

      // 4. Run quantification
      $('run-samples')?.click();
      setStatus('');
    } catch (err) {
      setStatus('');
      showErrors([`Test set loading failed: ${err.message || err}. Check network access to ENA and try again.`]);
    }
  }

  $('load-testset-yeast')?.addEventListener('click', () => loadTestSet('yeast-ampb'));
  $('load-testset-cneo')?.addEventListener('click', () => loadTestSet('cneo-lkh1'));

  init();

  return { loadTestSet, init };
}
