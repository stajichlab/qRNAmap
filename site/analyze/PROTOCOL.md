# Page worker protocol (`site/analyze/worker.js`)

The page has one module worker for all steps. `genome.js` starts it, and `app.js` gives it to the later steps as `window.qrnamapApp.getWorker()`.
The worker holds one thing after a genome is ready: the kallisto index (with the genome label and accession). The index never leaves the worker.
`tx2gene` and the gene table go to the page in `genome-ready`, and the page keeps them (`state.tx2gene`, `state.genesTsv`); the worker does not keep them.
The worker does not keep the transcript FASTA, the genome bytes or the GFF3 after the index is built. Task 6 needs only the index, `tx2gene` and the gene table (plan, Task 6).
kallisto runs through the vendored `ToolRunner` (`vendor/kallisto/runtime/index.mjs`). `ToolRunner` starts one more worker for each kallisto command and terminates it when the command ends.

## Rules

- Every message from the page has a `type` and a number `id`. Every answer has the `id` of its message.
- Take every id from `window.qrnamapApp.nextId()`. The genome step uses the same counter, so two messages never have the same id.
- The page keeps one pending request. It ignores an answer whose `id` is not the pending `id` (an answer to a cancelled run).
- Cancel: the page terminates the worker and starts a new one. The genome state of the old worker is lost. The page then sets `indexReady` to false.
  There is no cancel message, because a running WebAssembly call cannot stop by itself.
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
| `status` | `id` | Answers `status`. |

## Messages from the worker

| `type` | Fields | When |
|---|---|---|
| `progress` | `id`, `phase` (`read`, `build` or `index`), `text`, `done`, `total` (bytes; 0 when not known) | During a `genome` request. |
| `genome-ready` | `id`, `label`, `accession`, `nTranscripts`, `nGenes`, `bases` (sequence letters counted in the FASTA), `indexBytes`, `indexSeconds`, `indexArgs` (the kallisto argument list), `tx2gene` (text), `genesTsv` (text) | The index is built and kept in the worker. |
| `error` | `id`, `kind`, `message` | A `genome` request failed. The worker has no genome after an error. |
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

## Page state for the later steps

`window.qrnamapApp.state` (the same object as `state` in `genome.js`):

- `phase`: `choose`, `download`, `read`, `build`, `index` or `ready`. The section `#step-genome` has the same value in `data-phase`.
- `genome`: null, or `{ label, accession, source: 'catalog' | 'upload', bases }`.
- `nTranscripts`, `nGenes`: numbers (0 when no genome).
- `indexReady`: true only while the worker holds the index of `genome`.
- `tx2gene`, `genesTsv`: the texts from `genome-ready` (empty when no genome). Only the page keeps them.
- `indexArgs`: the kallisto index argument list (null when no genome).
- `downloadNotes`: catalog files whose downloaded size is not the catalog size (a note, not an error).

`window.qrnamapApp` = `{ state, getWorker, cancelGenome, nextId }`.

The page sends the event `qrnamap:genome` on `document` (with a copy of the state in `detail`) each time the state changes.

## For Task 6 (quantification)

Task 6 adds its own message types to this worker (for example `quant`), with ids from `nextId()`. The worker then writes the index from its genome state into the
`ToolRunner` input files. Gene counts need `tx2gene`: the page sends it with the request, or computes gene counts itself. A crash of the kallisto worker does not end the page worker, so the index stays. A crash of the page worker loses the index.
