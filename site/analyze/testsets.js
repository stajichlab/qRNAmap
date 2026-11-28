// Test sets: load yeast AmpB (GSE80357) and C. neoformans H99 Lkh1 (PRJNA976751).
import { VerifyingGunzip } from '../gz.js?v=__BUILD__';
import { firstNameStem } from './fastq.js?v=__BUILD__';

// The cut test-set reads (<run>_1.fastq.gz, <run>_2.fastq.gz) are GitHub release testsets-v1, mirrored to the lab S3 bucket
// with scripts/mirror-indexes-to-s3.sh testsets-v1. They are not in the repository or the site, so a clone and the desktop
// app stay small. A release download URL has no CORS header, so the page reads the S3 copy; when that fails it streams the
// run from ENA. `?testset_base=URL` replaces the S3 folder (other mirrors, a local copy).
const TESTSET_READS_URL = (new URLSearchParams(location.search).get('testset_base') || 'https://s3-west.nrp-nautilus.io/stajichlab/qRNAmap/testsets-v1').replace(/\/+$/, '');

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

export async function fetchTestsetRead(localUrl, enaUrl, maxRecords, signal) {
  if (localUrl) {
    try {
      const res = await fetch(localUrl, { signal });
      if (res.ok) {
        const buf = new Uint8Array(await res.arrayBuffer());
        const isGz = buf.length >= 2 && buf[0] === 0x1f && buf[1] === 0x8b;
        const chunks = [];
        let newlines = 0;
        const target = maxRecords * 4;
        let cut = false;

        if (isGz) {
          const gunzip = new VerifyingGunzip((chunk) => {
            if (cut) return;
            for (let i = 0; i < chunk.length; i++) {
              if (chunk[i] === 10) {
                newlines++;
                if (newlines === target) {
                  chunks.push(chunk.slice(0, i + 1));
                  cut = true;
                  return;
                }
              }
            }
            chunks.push(chunk);
          });
          gunzip.push(buf, true);
        } else {
          for (let i = 0; i < buf.length; i++) {
            if (buf[i] === 10) {
              newlines++;
              if (newlines === target) {
                chunks.push(buf.slice(0, i + 1));
                cut = true;
                break;
              }
            }
          }
          if (!cut) chunks.push(buf);
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
    } catch (_) {
      // Fall through to ENA stream
    }
  }
  return await streamCutFastq(enaUrl, maxRecords, signal);
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

      // 2. Immediately set up the table with all 6 sample rows
      const covName = ts.design?.replace(/^~\s*/, '') || 'condition';
      const initialRows = ts.samples.map((s) => {
        const covVal = s[covName] || s.condition || s.genotype || '';
        return {
          name: s.sample,
          values: [covVal],
          files: { 1: null, 2: null },
          expect: { 1: `${s.run}_1.fastq`, 2: `${s.run}_2.fastq` },
          loading: { 1: 'Waiting to download...', 2: 'Waiting to download...' },
          statusText: 'Waiting for download...',
        };
      });
      if (app.samplesStep?.setSamples) {
        app.samplesStep.setSamples({ columns: [covName], rows: initialRows });
      }

      // 3. Start catalog genome download and indexing concurrently
      setStatus(`Downloading and indexing catalog genome ${ts.catalog_accession} (${ts.organism})...`);
      const genomePromise = app.loadCatalogAccession ? app.loadCatalogAccession(ts.catalog_accession) : Promise.resolve();

      // 4. Download and process sample reads
      for (let i = 0; i < ts.samples.length; i++) {
        const s = ts.samples[i];
        app.samplesStep?.updateRow?.(i, {
          loading: { 1: 'Downloading mate 1...', 2: 'Downloading mate 2...' },
          statusText: 'Downloading reads...',
        });
        setStatus(`Fetching reads for sample ${i + 1} of ${ts.samples.length} (${s.sample}, ${s.run})...`);

        const local1 = `${TESTSET_READS_URL}/${s.run}_1.fastq.gz`;
        const local2 = `${TESTSET_READS_URL}/${s.run}_2.fastq.gz`;

        // Check for test mock override function
        const url1 = window.__testSetUrl ? window.__testSetUrl(s.run, 1) : s.fastq1_url;
        const url2 = window.__testSetUrl ? window.__testSetUrl(s.run, 2) : s.fastq2_url;

        const [m1, m2] = await Promise.all([
          fetchTestsetRead(local1, url1, ts.n_pairs),
          fetchTestsetRead(local2, url2, ts.n_pairs),
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

        app.samplesStep?.updateRow?.(i, {
          files: { 1: file1, 2: file2 },
          loading: null,
          statusText: 'Ready',
        });
      }

      // 5. Ensure genome indexing is complete
      setStatus('Waiting for genome index build to complete...');
      await genomePromise;

      setStatus('All samples ready. Counting reads...');

      // 6. Run quantification
      $('run-samples')?.click();
      setStatus('');
    } catch (err) {
      setStatus('');
      showErrors([`Test set loading failed: ${err.message || err}. Check network access and try again.`]);
    }
  }

  $('load-testset-yeast')?.addEventListener('click', () => loadTestSet('yeast-ampb'));
  $('load-testset-cneo')?.addEventListener('click', () => loadTestSet('cneo-lkh1'));

  init();

  return { loadTestSet, init };
}
