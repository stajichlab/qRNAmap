# Files copied from qmutfind

Source: https://github.com/stajichlab/qmutfind, commit `4404e1c93f47ea527287b93d6b496070563b193a`
(the head of the local qmutfind checkout on 2026-11-22). Both projects use GPL-3.0-only and have the same author.
Each copied file has a first line with the commit. That line is the only change unless the list below says more.

| File | Changes |
|---|---|
| `gz.js` | Header line only. |
| `reads.js` | Header line only. Not used by the page and not copied to `dist-site` (the build script omits it). It stays in `site/` because the plan names it. |
| `reference-panel.js` | Header line. The FungiDB command shown to the user is now `qrnamap fetch --source fungidb --id "<name>" --out genome/ --eupathdb-key-file ~/.eupathdb_key` (was `qmutfind fetch --organism ...`). The comment on line 187 names `qrnamap fetch`. The default genome list path (`genomes/genomes.json`) is not changed yet: a later task decides where the list comes from. |
| `persist.js` | Header line. The IndexedDB database name is `qrnamap` (was `qmutfind`), in the code and in the comment. |
| `theme.js` | Header line. The keys `qmutfind-theme` and `qmutfind-layout` and the event name `qmutfind-theme` are now `qrnamap-theme`, `qrnamap-layout` and `qrnamap-theme`. The page uses only the theme buttons (`light`, `dark`). The layout code stays, and does nothing without layout buttons. |
| `vendor/fflate.js` | Header line only. |
| `vendor/fflate.LICENSE` | No change (MIT, Arjun Barrett). |

`gz.js` and `reads.js` import `./vendor/fflate.js?v=__BUILD__` and `./gz.js?v=__BUILD__`. The build script replaces the token (as in qmutfind).
Measurements (docs/browser/page-v1-measurements.md) show that the gzip repack in `reads.js` does not meet the memory rule at the chosen preview size, so the page may not use `reads.js`. It is copied because the plan names it.

Open items:
- `reference-panel.js` still reads `genomes/genomes.json` by default. That file is not in `dist-site`. The genome step (`analyze/genome.js`) does not use the panel: it imports only `checkUrl` and `mb` from `reference-panel.js`, and its genome list comes from `catalog()` of qrnamap-wasm.
- The landing page links to `docs/browser/page-v1-measurements.md` on the `main` branch. That link returns 404 until branch `page-v1` is merged into `main`.
