#!/usr/bin/env python3
"""Merge the page-v1 measurement TSVs of one workflow run and apply the stated limit rules.

Usage: summarize.py RUN_DIR OUT_DIR
RUN_DIR holds the downloaded artifacts (any depth): quant.<engine>.<mode>.tsv, truncated.<engine>.tsv, genomes.<engine>.tsv.
Writes to OUT_DIR: quant_runs.tsv, quant_summary.tsv, limits.tsv, truncated.tsv, genomes_runs.tsv, genomes_summary.tsv.

Rules (from the Task 3 brief; the values are inputs, not results):
  time   per sample = prep_s + quant_s (the index is built once per genome, so it is not in this sum) <= 300 s
  memory N:          max(prep_tree_pss_peak_kb, quant_tree_pss_peak_kb) (whole browser) <= 1.5 GiB = 1,572,864 kB
  memory hard limit: the same peak <= 2 GiB = 2,097,152 kB
  a count passes for an engine only when all its runs finished and every run meets the rule.
  N and the hard limit are the largest passing count in all three browsers (Chromium, Firefox, WebKit), per mode.
"""
import csv, glob, os, sys

run_dir, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
BROWSERS = ['chromium', 'firefox', 'webkit']
TIME_S = 300.0
N_KB = 1572864
HARD_KB = 2097152


def read(pattern):
    rows = []
    for path in sorted(glob.glob(os.path.join(run_dir, '**', pattern), recursive=True)):
        with open(path, encoding='utf-8') as fh:
            rows.extend(csv.DictReader(fh, delimiter='\t'))
    return rows


def write(name, rows, cols=None):
    if not rows:
        open(os.path.join(out, name), 'w').close()
        return
    cols = cols or list(rows[0].keys())
    with open(os.path.join(out, name), 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, cols, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def fmt(values, digits=2):
    return ', '.join('' if v is None else (f'{v:.{digits}f}' if digits else str(int(v))) for v in values)


quant = read('quant.*.tsv')
write('quant_runs.tsv', quant)
groups = {}
for r in quant:
    groups.setdefault((r['engine'], r['mode'], int(r['pairs'])), []).append(r)
summary = []
for (engine, mode, pairs), rs in sorted(groups.items()):
    rs.sort(key=lambda r: int(r['rep']))
    sample_s = [None if num(r['prep_s']) is None or num(r['quant_s']) is None else num(r['prep_s']) + num(r['quant_s']) for r in rs]
    # Peak memory of one sample = the higher of the preparation peak and the quant peak (browsers); Node: maxRSS of the process.
    def peak_of(r):
        if engine == 'node':
            return num(r['node_maxrss_kb'])
        vals = [v for v in (num(r['prep_tree_pss_peak_kb']), num(r['quant_tree_pss_peak_kb'])) if v is not None]
        return max(vals) if vals else None
    peak = [peak_of(r) for r in rs]
    finished = sum(r['outcome'] == 'finished' for r in rs)
    ok_time = finished == len(rs) and all(v is not None and v <= TIME_S for v in sample_s)
    ok_n = ok_time and all(v is not None and v <= N_KB for v in peak)
    ok_hard = finished == len(rs) and all(v is not None and v <= HARD_KB for v in peak)
    summary.append({
        'engine': engine, 'mode': mode, 'pairs': pairs, 'runs': len(rs), 'finished': finished,
        'outcomes': ', '.join(r['outcome'] for r in rs),
        'index_s': fmt([num(r['index_s']) for r in rs]), 'prep_s': fmt([num(r['prep_s']) for r in rs]),
        'quant_s': fmt([num(r['quant_s']) for r in rs]), 'sample_s': fmt(sample_s),
        'peak_probe': 'node_maxrss_kb' if engine == 'node' else 'max(prep, quant) tree_pss_peak_kb',
        'prep_peak_mib': fmt([None if num(r['prep_tree_pss_peak_kb']) is None else num(r['prep_tree_pss_peak_kb']) / 1024 for r in rs], 1),
        'quant_peak_mib': fmt([None if num(r['quant_tree_pss_peak_kb']) is None else num(r['quant_tree_pss_peak_kb']) / 1024 for r in rs], 1),
        'peak_mib': fmt([None if v is None else v / 1024 for v in peak], 1),
        'baseline_mib': fmt([None if num(r['quant_tree_pss_baseline_kb']) is None else num(r['quant_tree_pss_baseline_kb']) / 1024 for r in rs], 1),
        'quant_heap_mib': fmt([None if num(r['quant_heap_final_bytes']) is None else num(r['quant_heap_final_bytes']) / 1048576 for r in rs], 1),
        'input_mib_per_mate': fmt([None if num(r['input1_bytes']) is None else num(r['input1_bytes']) / 1048576 for r in rs], 1),
        'n_pseudoaligned': fmt([num(r['n_pseudoaligned']) for r in rs], 0),
        'abundance_sha256_distinct': len({r['abundance_sha256'] for r in rs if r['abundance_sha256']}),
        'meets_time': ok_time, 'meets_n_rule': ok_n, 'meets_hard_rule': ok_hard,
        'details': ' | '.join(sorted({r['detail'][:200] for r in rs if r['detail']})),
    })
write('quant_summary.tsv', summary)

limits = []
by = {(s['engine'], s['mode'], s['pairs']): s for s in summary}
for mode in sorted({s['mode'] for s in summary}):
    counts = sorted({s['pairs'] for s in summary if s['mode'] == mode})
    for rule in ('meets_n_rule', 'meets_hard_rule'):
        passing = [c for c in counts if all(by.get((e, mode, c), {}).get(rule) for e in BROWSERS)]
        best = max(passing) if passing else ''
        first_fail = {e: next((c for c in counts if (e, mode, c) in by and not by[(e, mode, c)][rule]), '') for e in BROWSERS}
        limits.append({'mode': mode, 'rule': rule, 'largest_count_passing_all_browsers': best,
                       **{f'first_count_failing_{e}': first_fail[e] for e in BROWSERS}})
write('limits.tsv', limits)

# Same reads give the same abundance.tsv? One row per count: distinct SHA-256 over all engines and modes.
agree = []
for pairs in sorted({int(r['pairs']) for r in quant}):
    shas = {r['abundance_sha256'] for r in quant if int(r['pairs']) == pairs and r['outcome'] == 'finished'}
    agree.append({'pairs': pairs, 'finished_runs': sum(1 for r in quant if int(r['pairs']) == pairs and r['outcome'] == 'finished'),
                  'distinct_abundance_sha256': len(shas)})
write('quant_abundance_agreement.tsv', agree)

write('truncated.tsv', read('truncated.*.tsv'))
write('determinism.tsv', read('determinism.*.tsv'))

genomes = read('genomes.*.tsv')
write('genomes_runs.tsv', genomes)
gg = {}
for r in genomes:
    gg.setdefault((r['accession'], r['engine']), []).append(r)
gsum = []
for (acc, engine), rs in sorted(gg.items()):
    rs.sort(key=lambda r: int(r['rep']))
    first = rs[0]
    gsum.append({
        'accession': acc, 'engine': engine, 'runs': len(rs), 'finished': sum(r['outcome'] == 'finished' for r in rs),
        'fna_gz_bytes': first['fna_gz_bytes'], 'gff_gz_bytes': first['gff_gz_bytes'], 'fna_bytes': first['fna_bytes'],
        'n_transcripts': first['n_transcripts'], 'n_genes': first['n_genes'], 'transcripts_bytes': first['transcripts_bytes'],
        'transcripts_sha256_distinct': len({r['transcripts_sha256'] for r in rs if r['transcripts_sha256']}),
        'gunzip_fflate_s': fmt([num(r['gunzip_fflate_s']) for r in rs]), 'gunzip_ds_s': fmt([num(r['gunzip_ds_s']) for r in rs]),
        'gunzip_equal': ', '.join(r['gunzip_equal'] for r in rs), 'members': f"{first['fna_members']}/{first['gff_members']}",
        'build_s': fmt([num(r['build_s']) for r in rs]), 'index_s': fmt([num(r['index_s']) for r in rs]),
        'index_heap_mib': fmt([None if num(r['index_heap_final_bytes']) is None else num(r['index_heap_final_bytes']) / 1048576 for r in rs], 1),
        'build_peak_mib': fmt([None if num(r['build_tree_pss_peak_kb'] or r['node_rss_build_peak_kb']) is None else num(r['build_tree_pss_peak_kb'] or r['node_rss_build_peak_kb']) / 1024 for r in rs], 1),
        'index_peak_mib': fmt([None if num(r['index_tree_pss_peak_kb'] or r['node_rss_index_peak_kb']) is None else num(r['index_tree_pss_peak_kb'] or r['node_rss_index_peak_kb']) / 1024 for r in rs], 1),
        'node_maxrss_mib': fmt([None if num(r['node_maxrss_kb']) is None else num(r['node_maxrss_kb']) / 1024 for r in rs], 1),
        'ncbi_fetch_s': first['ncbi_fetch_s'], 'ncbi_fetch_bytes': first['ncbi_fetch_bytes'], 'ncbi_sha_equal': first['ncbi_sha_equal'],
        'details': ' | '.join(sorted({r['detail'][:200] for r in rs if r['detail']})),
    })
write('genomes_summary.tsv', gsum)
print(f'quant rows {len(quant)}, genome rows {len(genomes)}')
