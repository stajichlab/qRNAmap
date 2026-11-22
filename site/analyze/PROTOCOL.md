# Page worker protocol (`site/analyze/worker.js`)

The page has one module worker for all steps. `genome.js` starts it, and `app.js` gives it to the later steps as `window.qrnamapApp.getWorker()`.
`genome.js` handles the `genome` request; `quant.js` handles the `quant` requests (it adds its own `message` listener and ignores other ids).
The worker holds one thing after a genome is ready: the kallisto index (with the genome label and accession). The index never leaves the worker.
`tx2gene` and the gene table go to the page in `genome-ready`, and the page keeps them (`state.tx2gene`, `state.genesTsv`); the worker does not keep them.
The worker does not keep the transcript FASTA, the genome bytes or the GFF3 after the index is built. Task 6 needs only the index, `tx2gene` and the gene table (plan, Task 6).
kallisto runs through the vendored `ToolRunner` (`vendor/kallisto/runtime/index.mjs`). `ToolRunner` starts one more worker for each kallisto command and terminates it when the command ends.

## Rules

- Every message from the page has a `type` and a number `id`. Every answer has the `id` of its message.
- Take every id from `window.qrnamapApp.nextId()`. The genome step uses the same counter, so two messages never have the same id.
- The page keeps one pending request. It ignores an answer whose `id` is not the pending `id` (an answer to a cancelled run).
- Cancel of the genome step: the page terminates the worker and starts a new one. The genome state of the old worker is lost. The page then sets `indexReady` to false.
  There is no cancel message for a `genome` request, because a running WebAssembly call in the worker cannot stop by itself.
- Cancel of a `quant` request: the page sends `quant-cancel` with the id of the request. The worker stops reading the files, or aborts the
  `ToolRunner` run (which terminates the kallisto worker), and answers `error` with kind `cancelled`. The page worker and its index stay,
  so a new run needs no new index. kallisto runs in its own worker, so the page worker can read the cancel message while kallisto runs.
- "Choose another genome" also replaces the worker, so the old index is dropped.
- Each time the worker is replaced, the page sends the event `qrnamap:worker-replaced` on `window` (a `CustomEvent`, `detail: { reason }`,
  `reason` is `cancel`, `crash` or `change`). A step that has a request open on the old worker must treat it as lost, and must take the new worker from `getWorker()`.
- A crash of the worker (the `error` or `messageerror` event, for example when the browser runs out of memory) rejects the pending request with
  kind `memory`. The page then starts a new worker. If no request was pending and an index was ready, the page says that the index is lost.
- A new `genome` message clears the old genome at once, before any work. So no step can use the index of the previous genome.
- Memory: the worker sets the request fields `fasta` and `gff` to null when it starts to read them, and reads and builds in a helper, so the genome bytes,
  the GFF3 text and the File objects are unreachable while kallisto runs. The page also drops its references to the files before it waits for the answer.
  The vendored `ToolRunner` keeps the transcript FASTA until kallisto returns (it copies it to the kallisto worker). Peak memory of this path: UNVERIFIED
  (docs/browser/page-v1-measurements.md, "Page genome path: peak memory").

## Messages from the page

| `type` | Fields | What the worker does |
|---|---|---|
| `genome` | `id`, `fasta` (File or Blob), `gff` (File or Blob), `label`, `accession` (empty for uploads), `limitBases` (number), `limitText` (text added to the size message) | Reads both files whole. A file that starts with the gzip bytes `1f 8b` is unpacked through `gz.js` (fflate, with the CRC-32 and length check of every member). Counts the sequence letters of the FASTA and stops when the count is above `limitBases`. Calls `build_transcriptome` (qrnamap-wasm). Runs `kallisto index` with the argument list from `kallisto_args({"kind":"index","fasta":"transcripts.fa","index":"transcripts.idx"})`, which is `["index","-i","transcripts.idx","transcripts.fa"]` (the same list as the Task 3 measurements; no `-t`). |
| `quant` | `id`, `name` (sample name), `mate1` (File), `mate2` (File or null for single-end), `pairs` (N, 1 to 200,000), `fragment` (`[mean, sd]` for single-end, else null), `nTranscripts` (from `genome-ready`) | Needs the index (else `error` kind `no-index`). Reads the first `pairs` records of each mate (`fastq.js`; cut rule below), checks that the kept bytes are FASTQ, checks that both mates have the same record count, runs `kallisto quant` with the argument list from `kallisto_args({"kind":"quant","sample":{...},"index":"transcripts.idx","out":"kq","threads":1,...})` and the files `transcripts.idx`, `r1.fastq` (and `r2.fastq`), checks that `n_processed` equals the pairs given, and checks `abundance.tsv` (`abundance-check.js`). When that check fails, it runs kallisto once more with the same bytes; a second failure is `error` kind `output`. Then it answers `quant-done`. The requests run one at a time, in the order they came. Nothing of a sample is kept after the answer. |
| `quant-cancel` | `id` (the id of the `quant` request) | Stops that request (see "Cancel of a `quant` request"). |
| `status` | `id` | Answers `status`. |

## Messages from the worker

| `type` | Fields | When |
|---|---|---|
| `progress` | `id`, `phase` (`read`, `build` or `index`), `text`, `done`, `total` (bytes; 0 when not known) | During a `genome` request. |
| `progress` | `id`, `phase` (`prep` or `quant`), `text` | During a `quant` request. |
| `quant-done` | `id`, `name`, `nPairs` (records per mate given to kallisto), `wholeFile` (true when the files had fewer than `pairs` records: every record was used), `paired`, `args` (the kallisto argument list), `reruns` (0 or 1), `rerunReason` (the problems of the first output, or null), `abundance` (the text of `abundance.tsv`), `runInfo` (`{n_processed, n_pseudoaligned, p_pseudoaligned}`) | The sample is counted. |
| `genome-ready` | `id`, `label`, `accession`, `nTranscripts`, `nGenes`, `bases` (sequence letters counted in the FASTA), `indexBytes`, `indexSeconds`, `indexArgs` (the kallisto argument list), `tx2gene` (text), `genesTsv` (text) | The index is built and kept in the worker. |
| `error` | `id`, `kind`, `message` | A `genome` request failed (the worker has no genome after it), or a `quant` request failed (the index stays). |
| `status` | `id`, `genome` (null, or `{label, accession, indexBytes}`) | Answer to `status`. |

Error kinds:

| `kind` | Cause | What the page adds |
|---|---|---|
| `input` | An empty file, or a `build_transcriptome` error (for example FASTA names that do not match the GFF3, or a GFF3 with no usable transcripts). The message holds the text of `build_transcriptome`. | Nothing. |
| `gzip` | The gzip data is cut short or damaged (fflate or `gz.js` message, for example `unexpected EOF`). | Nothing. |
| `limit` | The FASTA has more sequence letters than `limitBases`, or the GFF3 unpacks to more than 260 MB (a guard: about 4 times the largest GFF3 of the Task 3 runs, 64,659,993 bytes). The FASTA size message ends with `limitText`, which names the command-line route. | Nothing (the route is in the message). |
| `kallisto` | `kallisto index` exited with a code other than 0, or wrote no index. | The command-line route. |
| `memory` | `ToolRunner` failed with an error that looks like a memory failure (`RangeError`, "out of memory", "could not allocate", and similar), or the page worker crashed. | The command-line route. |
| `internal` | Any other error. | Nothing. |
| `no-index` | (`quant`) The worker has no index. | Nothing. |
| `cancelled` | (`quant`) The page sent `quant-cancel`. | Nothing. |
| `output` | (`quant`) The `abundance.tsv` failed the check twice, or `run_info.json` cannot be read. | Nothing. |

For a `quant` request, `input` also covers the read files: an empty file, a gzip file that is cut short or damaged, a file that ends inside a record,
a record that is not FASTQ (no `@` line, no `+` line, or a quality line with another length than the sequence), and mates with different record counts.
`kallisto` also covers `n_processed` different from the pairs given.

## Read cut rule (Task 9 needs it)

- For each mate: stream the File (`File.stream()`). A file that starts with the bytes `1f 8b` goes through `gz.js` (fflate, every gzip member, CRC-32 and length check);
  any other file is read as plain FASTQ. Count the newline bytes (10) of the plain bytes. Keep the bytes up to and with newline number `pairs x 4`. Cancel the stream.
- If the file ends first: the gzip end is checked (a cut gzip file is an error); a last line without a newline gets one; the line count must be a multiple of 4.
  Then every record of the file is used, and `wholeFile` is true.
- Both mates must have the same record count. When both files have more than `pairs` records, both are cut at `pairs`, so a difference after record `pairs`
  is not seen (the page does not read the rest of the files).
- kallisto gets the kept bytes as plain FASTQ (`r1.fastq`, `r2.fastq`), never gzip. The page does not refuse a file by its size: it reads only the first `pairs` records.
- `DecompressionStream('gzip')` is not used: in Chromium 156 it fails with "Junk found after end of compressed data" on a gzip file with two members (2026-11-22).
- Memory: the kept bytes are held in the page worker, copied once to the kallisto worker by `ToolRunner` (structured clone), and once more into MEMFS
  (about 3 copies, as in the measurements). The index is also copied to each kallisto run. Peak memory of this page path: UNVERIFIED (not measured).

## Page state for the later steps

`window.qrnamapApp.state` (the same object as `state` in `genome.js`):

- `phase`: `choose`, `download`, `read`, `build`, `index` or `ready`. The section `#step-genome` has the same value in `data-phase`.
- `genome`: null, or `{ label, accession, source: 'catalog' | 'upload', bases }`.
- `nTranscripts`, `nGenes`: numbers (0 when no genome).
- `indexReady`: true only while the worker holds the index of `genome`.
- `tx2gene`, `genesTsv`: the texts from `genome-ready` (empty when no genome). Only the page keeps them.
- `indexArgs`: the kallisto index argument list (null when no genome).
- `downloadNotes`: catalog files whose downloaded size is not the catalog size (a note, not an error).
- `samples` (step 2): one object per sample of the last run, in table order: `{ name, covariates ({column: value}, as in the sheet; a numeric column has no
  `:numeric` here), status ('waiting', 'running', 'done', 'error', 'cancelled' or 'not run'), nPairs (records per mate given to kallisto), processed
  (n_processed), pseudoaligned (n_pseudoaligned), fraction (pseudoaligned / nPairs), allPairs (true when the files had fewer than N pairs), geneCounts
  (object {gene: integer}, from `gene_counts` of qrnamap-wasm with `tx2gene`; null until done), error (text), kind (error kind), reruns (0 or 1), args (the kallisto argument list) }`.
  Empty when no run was made, and cleared when the worker is replaced (the counts belong to the old index).
- `preview` (step 2): `{ n, label }`, for example `{ n: 200000, label: 'preview of the first 200,000 pairs' }`. Every result made from the samples must carry
  `label`; a sample with `allPairs` true says "all pairs (X)" instead.
- `samplesPhase` (step 2): 'edit', 'running', 'done', 'cancelled' or 'lost'. `#step-samples` has the same value in `data-phase` (`edit` after a clear).

`window.qrnamapApp` = `{ state, getWorker, cancelGenome, nextId }`.

The page sends the event `qrnamap:genome` on `document` (with a copy of the state in `detail`) each time the state changes.
The page sends the event `qrnamap:samples-done` on `document` when a run of all samples has ended without a cancel (`detail: { samples, preview }`).
Samples with an error are in the list with status `error`. Step 2 is marked done when at least one sample is `done`.

## Steps 1 and 2

- The samples step opens when the index is ready (`completeStep('genome')`) and locks again when the genome is cleared (`resetFrom('genome')`).
- After a run with at least one counted sample, step 2 is marked done with `markDone('samples')`. Step 3 stays locked until Task 7 builds it.
