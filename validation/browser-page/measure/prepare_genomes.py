#!/usr/bin/env python3
"""Page-v1 measurements, CI only: download the genomes to measure from NCBI (URLs from the catalog) and time the downloads.

Usage: prepare_genomes.py CATALOG_JSON OUT_DIR
Genomes: S. cerevisiae S288C, A. fumigatus Af293, C. neoformans H99 and the three largest genomes (genome_size) of the
catalog that have a GFF3. Writes OUT_DIR/genomes.tsv (the measure.mjs --genomes table) and OUT_DIR/download.tsv
(bytes, catalog bytes, seconds, SHA-256 per file). Stdlib only.
"""
import hashlib, json, os, sys, time, urllib.request

FIXED = ['GCF_000146045.2', 'GCF_000002655.1', 'GCF_000149245.1']

catalog, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
genomes = json.load(open(catalog, encoding='utf-8'))['genomes']
usable = [g for g in genomes if g.get('gff')]
largest = [g['accession'] for g in sorted(usable, key=lambda g: -g['genome_size'])[:3]]
by_acc = {g['accession']: g for g in usable}
chosen = FIXED + [a for a in largest if a not in FIXED]

rows, dl = [], []
for acc in chosen:
    g = by_acc[acc]
    names = {}
    for kind in ('fasta', 'gff'):
        ref = g[kind]
        path = os.path.join(out, ref['name'])
        for attempt in range(3):
            try:
                start = time.monotonic()
                req = urllib.request.Request(ref['url'], headers={'User-Agent': 'qrnamap-measure'})
                with urllib.request.urlopen(req, timeout=600) as resp, open(path, 'wb') as fh:
                    while True:
                        chunk = resp.read(1 << 20)
                        if not chunk:
                            break
                        fh.write(chunk)
                secs = time.monotonic() - start
                break
            except Exception as exc:  # retry, then fail
                if attempt == 2:
                    raise
                print(f'retry {ref["url"]}: {exc}', file=sys.stderr)
                time.sleep(10)
        data = open(path, 'rb').read()
        sha = hashlib.sha256(data).hexdigest()
        names[kind] = (ref['name'], ref['url'], sha)
        dl.append([acc, g['label'], kind, ref['name'], str(len(data)), str(ref['bytes']), f'{secs:.2f}', sha])
    rows.append([acc, names['fasta'][0], names['gff'][0], names['fasta'][1], names['gff'][1], names['fasta'][2], names['gff'][2],
                 g['label'], str(g['genome_size']), str(g['genes']), 'largest3' if acc in largest else 'named'])

with open(os.path.join(out, 'genomes.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('accession\tfna\tgff\tfna_url\tgff_url\tfna_sha256\tgff_sha256\tlabel\tgenome_size\tgenes\twhy\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
with open(os.path.join(out, 'download.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('accession\tlabel\tkind\tfile\tbytes\tcatalog_bytes\tdownload_s\tsha256\n')
    for r in dl:
        fh.write('\t'.join(r) + '\n')
print(open(os.path.join(out, 'download.tsv'), encoding='utf-8').read())
