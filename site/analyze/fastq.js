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
 * One pass over the first `maxRecords` records of `file`. `sink(bytes)` gets every kept chunk (plain bytes) in order.
 * Returns `{ total, lines, cut, addNewline }`. Checks the gzip end (when the whole file was read) and blank lines.
 */
async function scan(file, maxRecords, label, signal, onRecords, sink, onChunk) {
  const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  const gz = head[0] === 0x1f && head[1] === 0x8b;
  const target = maxRecords * 4;
  let total = 0;
  let lines = 0;
  let cut = false;
  let reported = 0;
  let lineLength = 0; // bytes of the current line, without a carriage return
  let blank = 0; // the number of the first blank line (0: none)
  let last = -1; // the last kept byte
  const feed = (bytes) => {
    if (cut) return;
    let seen = 0;
    for (let i = 0; i < bytes.length; i++) {
      const b = bytes[i];
      if (b === NL) {
        if (lineLength === 0 && !blank) blank = lines + seen + 1;
        lineLength = 0;
        if (++seen + lines === target) { bytes = bytes.subarray(0, i + 1); cut = true; break; }
      } else if (b !== 13) lineLength++;
    }
    lines += seen;
    if (bytes.length) { last = bytes[bytes.length - 1]; sink(bytes); }
    total += bytes.length;
  };
  const gunzip = gz ? new VerifyingGunzip(feed) : null;
  const reader = file.stream().getReader();
  let lastYield = performance.now();
  try {
    while (!cut) {
      // Give the event loop a turn at least every 50 ms. Reason: when the file data is already in memory, an engine can resolve
      // every reader.read() as a microtask, and then a `quant-cancel` message (a task) would wait until the whole read ends.
      if (performance.now() - lastYield > 50) {
        await new Promise((done) => setTimeout(done, 0));
        lastYield = performance.now();
      }
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
      if (onChunk) await onChunk();
      if (lines - reported >= 40000) { reported = lines; onRecords(Math.floor(lines / 4)); }
    }
  } finally {
    reader.cancel().catch(() => {});
  }
  // A file whose last line has no newline: the line still counts.
  const addNewline = !cut && total > 0 && last !== NL;
  if (addNewline) lines += 1;
  if (total === 0) throw new ReadError(`${label} is empty after gunzip.`);
  if (blank) throw new ReadError(`${label} has a blank line (line ${blank.toLocaleString('en-US')}). A FASTQ file has no blank lines. Remove it, or check that the file is FASTQ.`);
  if (lines % 4 !== 0) {
    throw new ReadError(`${label} ends inside a record: it has ${lines} lines, which is not a multiple of 4. A FASTQ record has 4 lines. The file may be cut short.`);
  }
  return { total, lines, cut, addNewline };
}

/**
 * Keep the first `maxRecords` records (4 lines each) of `file` (a File or Blob). `label` names the file in messages.
 * Returns `{ bytes, records, wholeFile }`: `wholeFile` is true when the file ended before `maxRecords` records
 * (then every record of the file is kept). `signal` (AbortSignal) stops the read. `onRecords(n)` gets the records kept so far.
 * Memory: the file is read twice. Pass 1 only counts the bytes to keep. Pass 2 copies them into one array of that exact size.
 * So the kept bytes exist once (plus one chunk), not as a list of chunks and a joined copy at the same time. Cost: the reading
 * and gunzip time is about twice that of one pass (reasoning, not a measurement).
 */
export async function firstRecords(file, maxRecords, label, { signal, onRecords = () => {} } = {}) {
  if (file.size === 0) throw new ReadError(`${label} is empty (0 bytes).`);
  const first = await scan(file, maxRecords, label, signal, onRecords, () => {});
  const bytes = new Uint8Array(first.total + (first.addNewline ? 1 : 0));
  let offset = 0;
  const second = await scan(file, maxRecords, label, signal, () => {}, (chunk) => {
    if (offset + chunk.length > first.total) throw new ReadError(`${label} changed while it was read. Pick it again.`);
    bytes.set(chunk, offset);
    offset += chunk.length;
  });
  if (second.total !== first.total || second.lines !== first.lines) throw new ReadError(`${label} changed while it was read. Pick it again.`);
  if (first.addNewline) bytes[first.total] = NL;
  onRecords(first.lines / 4);
  return { bytes, records: first.lines / 4, wholeFile: !first.cut };
}

/**
 * Stream records of `file` directly into a `StreamRingBuffer`.
 * Decompresses and writes chunks on-the-fly, pausing via ring buffer backpressure.
 * Closes the ring buffer on completion.
 */
export async function streamFastqRecords(file, maxRecords, label, ring, { signal, onRecords = () => {} } = {}) {
  if (file.size === 0) {
    ring.abort();
    throw new ReadError(`${label} is empty (0 bytes).`);
  }
  const pending = [];
  const sink = (chunk) => {
    pending.push(chunk);
  };
  const onChunk = async () => {
    while (pending.length) {
      const chunk = pending.shift();
      await ring.write(chunk);
    }
  };
  try {
    const res = await scan(file, maxRecords, label, signal, onRecords, sink, onChunk);
    await onChunk();
    if (res.addNewline) {
      await ring.write(new Uint8Array([NL]));
    }
    ring.close();
    onRecords(Math.floor(res.lines / 4));
    return { records: Math.floor(res.lines / 4), wholeFile: !res.cut };
  } catch (error) {
    ring.abort();
    throw error;
  }
}


/** The sequence length of the first record of FASTQ bytes (line 2, without a carriage return). */
export function firstReadLength(bytes) {
  let line = 0;
  let start = 0;
  for (let i = 0; i < bytes.length; i++) {
    if (bytes[i] !== NL) continue;
    if (line === 1) return (i > start && bytes[i - 1] === 13) ? i - start - 1 : i - start;
    line++;
    start = i + 1;
  }
  return 0;
}

/** The first word of the first line of FASTQ bytes, without `@` and without a `/1` or `/2` end (for the mate name check). */
export function firstNameStem(bytes) {
  let end = 0;
  while (end < bytes.length && end < 4096 && bytes[end] !== NL && bytes[end] !== 13 && bytes[end] !== 32 && bytes[end] !== 9) end++;
  return new TextDecoder().decode(bytes.subarray(1, end)).replace(/\/[12]$/, '');
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
