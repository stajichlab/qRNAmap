# Page worker protocol (`site/analyze/worker.js`)

The page has one module worker for all steps. `genome.js` starts it, and `app.js` gives it to the later steps as `window.qrnamapApp.getWorker()`.
The worker holds the genome state: the transcript FASTA, `tx2gene`, the gene table and the kallisto index. The index never leaves the worker.
kallisto runs through the vendored `ToolRunner` (`vendor/kallisto/runtime/index.mjs`). `ToolRunner` starts one more worker for each kallisto command and terminates it when the command ends.

## Rules

- Every message from the page has a `type` and a number `id`. Every answer has the `id` of its message.
- The page keeps one pending request. It ignores an answer whose `id` is not the pending `id` (an answer to a cancelled run).
- Cancel: the page terminates the worker and starts a new one. The genome state of the old worker is lost. The page then sets `indexReady` to false.
  There is no cancel message, because a running WebAssembly call cannot stop by itself.
- A crash of the worker (the `error` or `messageerror` event, for example when the browser runs out of memory) rejects the pending request with
  kind `memory`. The page then starts a new worker. If no request was pending and an index was ready, the page says that the index is lost.
- A new `genome` message clears the old genome at once, before any work. So no step can use the index of the previous genome.

## Messages from the page

| `type` | Fields | What the worker does |
|---|---|---|
| `genome` | `id`, `fasta` (File or Blob), `gff` (File or Blob), `label`, `accession` (empty for uploads), `limitBases` (number), `limitText` (text added to the size message) | Reads both files whole. A file that starts with the gzip bytes `1f 8b` is unpacked through `gz.js` (fflate, with the CRC-32 and length check of every member). Counts the sequence letters of the FASTA and stops when the count is above `limitBases`. Calls `build_transcriptome` (qrnamap-wasm). Runs `kallisto index` with the argument list from `kallisto_args({"kind":"index","fasta":"transcripts.fa","index":"transcripts.idx"})`, which is `["index","-i","transcripts.idx","transcripts.fa"]` (the same list as the Task 3 measurements; no `-t`). |
| `clear` | `id` | Drops the genome state. No answer. |
| `status` | `id` | Answers `status`. |

## Messages from the worker

| `type` | Fields | When |
|---|---|---|
| `progress` | `id`, `phase` (`read`, `build` or `index`), `text`, `done`, `total` (bytes; 0 when not known) | During a `genome` request. |
| `genome-ready` | `id`, `label`, `accession`, `nTranscripts`, `nGenes`, `bases` (sequence letters counted in the FASTA), `indexBytes`, `indexSeconds`, `tx2gene` (text), `genesTsv` (text) | The index is built and kept in the worker. |
| `error` | `id`, `kind`, `message` | A `genome` request failed. The worker has no genome after an error. |
| `status` | `id`, `genome` (null, or `{label, accession, nTranscripts, nGenes, indexBytes}`) | Answer to `status`. |

Error kinds:

| `kind` | Cause | What the page adds |
|---|---|---|
| `input` | An empty file, or a `build_transcriptome` error (for example FASTA names that do not match the GFF3, or a GFF3 with no usable transcripts). The message holds the text of `build_transcriptome`. | Nothing. |
| `gzip` | The gzip data is cut short or damaged (fflate or `gz.js` message, for example `unexpected EOF`). | Nothing. |
| `limit` | The FASTA has more sequence letters than `limitBases`, or the GFF3 unpacks to more than 1 GB. The size message ends with `limitText`, which names the command-line route. | Nothing (the route is in the message). |
| `kallisto` | `kallisto index` exited with a code other than 0, or wrote no index. | The command-line route. |
| `memory` | `ToolRunner` failed with an error that looks like a memory failure (`RangeError`, "out of memory", "could not allocate", and similar), or the page worker crashed. | The command-line route. |
| `internal` | Any other error. | Nothing. |

## Page state for the later steps

`window.qrnamapApp.state` (the same object as `state` in `genome.js`):

- `phase`: `choose`, `download`, `read`, `build`, `index` or `ready`. The section `#step-genome` has the same value in `data-phase`.
- `genome`: null, or `{ label, accession, source: 'catalog' | 'upload', bases }`.
- `nTranscripts`, `nGenes`: numbers (0 when no genome).
- `indexReady`: true only while the worker holds the index of `genome`.
- `tx2gene`, `genesTsv`: the texts from `genome-ready` (empty when no genome).

The page sends the event `qrnamap:genome` on `document` (with a copy of the state in `detail`) each time the state changes.

## For Task 6 (quantification)

Task 6 adds its own message types to this worker (for example `quant`). The worker then writes the index from its genome state into the
`ToolRunner` input files. A crash of the kallisto worker does not end the page worker, so the index stays. A crash of the page worker loses the index.
