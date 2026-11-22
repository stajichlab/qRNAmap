// Read the first N records of a local FASTQ file (plain or gzip) for the preview, and check that the kept bytes are FASTQ.
// Cut rule (docs/browser/page-v1-measurements.md, "Truncated gzip from ENA"): stream the file, count newlines up to N x 4,
// keep the bytes up to that newline, cancel the stream. The rest of the file is never read.
// The gzip reader is gz.js (fflate plus the CRC-32 and length check of every member), as in the genome step. The built-in
// DecompressionStream is not used: in Chromium 156 it stops with "Junk found after end of compressed data" on a gzip file
// with two members (checked on 2026-11-22), and FASTQ files made with bgzip have many members.
// The decoder end is called only when the whole file was read (a cut stream would throw "unexpected EOF").
import { VerifyingGunzip } from '../gz.js?v=__BUILD__';

/** A problem with the input file. The message is for the user. */
export class ReadError extends Error {}

const NL = 10;

/**
 * Keep the first `maxRecords` records (4 lines each) of `file` (a File or Blob). `label` names the file in messages.
 * Returns `{ bytes, records, wholeFile }`: `wholeFile` is true when the file ended before `maxRecords` records
 * (then every record of the file is kept). `signal` (AbortSignal) stops the read. `onRecords(n)` gets the records kept so far.
 */
export async function firstRecords(file, maxRecords, label, { signal, onRecords = () => {} } = {}) {
  if (file.size === 0) throw new ReadError(`${label} is empty (0 bytes).`);
  const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  const gz = head[0] === 0x1f && head[1] === 0x8b;
  const target = maxRecords * 4;
  const parts = [];
  let total = 0;
  let lines = 0;
  let cut = false;
  let reported = 0;
  const feed = (bytes) => {
    if (cut) return;
    let seen = 0;
    for (let i = 0; i < bytes.length; i++) {
      if (bytes[i] === NL && ++seen + lines === target) { bytes = bytes.subarray(0, i + 1); cut = true; break; }
    }
    lines += seen;
    parts.push(bytes);
    total += bytes.length;
  };
  const gunzip = gz ? new VerifyingGunzip(feed) : null;
  const reader = file.stream().getReader();
  try {
    while (!cut) {
      if (signal?.aborted) throw signal.reason ?? new DOMException('cancelled', 'AbortError');
      const { done, value } = await reader.read();
      if (done) {
        if (gunzip) {
          try { gunzip.push(new Uint8Array(0), true); } catch (error) {
            throw new ReadError(`${label} is not a complete gzip file (${error.message}). The file is damaged or cut short.`);
          }
        }
        break;
      }
      if (gunzip) {
        try { gunzip.push(value, false); } catch (error) {
          if (!cut) throw new ReadError(`${label} is not a valid gzip file (${error.message}).`);
        }
      } else {
        feed(value);
      }
      if (lines - reported >= 40000) { reported = lines; onRecords(Math.floor(lines / 4)); }
    }
  } finally {
    reader.cancel().catch(() => {});
  }
  // A file whose last line has no newline: the line still counts.
  const last = parts.length ? parts[parts.length - 1] : null;
  if (!cut && last && last.length && last[last.length - 1] !== NL) {
    parts.push(new Uint8Array([NL]));
    total += 1;
    lines += 1;
  }
  if (total === 0) throw new ReadError(`${label} is empty after gunzip.`);
  if (lines % 4 !== 0) {
    throw new ReadError(`${label} ends inside a record: it has ${lines} lines, which is not a multiple of 4. A FASTQ record has 4 lines. The file may be cut short.`);
  }
  const bytes = new Uint8Array(total);
  let offset = 0;
  for (const p of parts) { bytes.set(p, offset); offset += p.length; }
  parts.length = 0;
  onRecords(lines / 4);
  return { bytes, records: lines / 4, wholeFile: !cut };
}

/**
 * Check that `bytes` (whole records, from firstRecords) are FASTQ: line 1 of each record starts with `@`, line 3 starts
 * with `+`, and the sequence and the quality have the same length. Throws a ReadError that names the first bad record.
 */
export function checkFastq(bytes, label) {
  let start = 0;
  let line = 0;
  let seqLength = 0;
  for (let i = 0; i < bytes.length; i++) {
    if (bytes[i] !== NL) continue;
    let end = i;
    if (end > start && bytes[end - 1] === 13) end--; // CRLF
    const kind = line % 4;
    const record = Math.floor(line / 4) + 1;
    const bad = (why) => new ReadError(`${label} is not a FASTQ file: record ${record.toLocaleString('en-US')} (line ${(line + 1).toLocaleString('en-US')}) ${why}.`);
    if (kind === 0 && bytes[start] !== 64) throw bad('does not start with "@"');
    if (kind === 1) seqLength = end - start;
    if (kind === 2 && bytes[start] !== 43) throw bad('has no "+" line');
    if (kind === 3 && end - start !== seqLength) throw bad(`has ${seqLength} bases but ${end - start} quality letters`);
    line++;
    start = i + 1;
  }
}
