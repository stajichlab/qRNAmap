// The WebAssembly package must give the same files as the native CLI for the same input.
// Run with `pixi run test-wasm` (builds pkg-node and the CLI first). No BASH_SOURCE is used.
//
// Inputs (the yeast files the browser spike 0 made, /bigdata, not in git):
//   genomes/yeast/GCF_000146045.2_R64_genomic.{fna,gff}.gz  (NCBI RefSeq, S. cerevisiae)
//   native/kallisto/abundance.tsv                            (native kallisto 0.52.0, one sample)
// Set QRNAMAP_FIXTURES to another directory with the same layout. The test is skipped when the
// directory is missing, unless QRNAMAP_REQUIRE_TOOLS=1.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, readFileSync, existsSync, rmSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const require = createRequire(import.meta.url);
const wasm = require(join(process.cwd(), 'pkg-node', 'qrnamap_wasm.js'));
const CLI = process.env.QRNAMAP_BIN;
const FIX = process.env.QRNAMAP_FIXTURES || '/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/fixtures';
const FASTA_GZ = join(FIX, 'genomes/yeast/GCF_000146045.2_R64_genomic.fna.gz');
const GFF_GZ = join(FIX, 'genomes/yeast/GCF_000146045.2_R64_genomic.gff.gz');
const ABUNDANCE = join(FIX, 'native/kallisto/abundance.tsv');
const have = [FASTA_GZ, GFF_GZ, ABUNDANCE].every(existsSync);
const skip = !have && process.env.QRNAMAP_REQUIRE_TOOLS !== '1' ? 'fixtures not found' : false;
if (!have && !skip) throw new Error(`fixtures missing under ${FIX}`);

function run(cmd, args) {
  const r = spawnSync(cmd, args, { encoding: 'buffer', maxBuffer: 1 << 30 });
  assert.equal(r.status, 0, `${cmd} ${args.join(' ')}: ${r.stderr}`);
  return r.stdout;
}

test('the catalog lists genomes', () => {
  const g = JSON.parse(wasm.catalog());
  assert.ok(g.length > 0);
  assert.ok(g.every((x) => x.accession && x.fasta.url && x.gff.url));
});

test('build_transcriptome gives the bytes of `qrnamap transcripts`', { skip }, () => {
  const dir = mkdtempSync(join(tmpdir(), 'qrnamap-wasm-'));
  try {
    run(CLI, ['transcripts', '--reference', FASTA_GZ, '--gff', GFF_GZ,
      '--out', join(dir, 't.fa'), '--tx2gene', join(dir, 't2g.tsv'), '--genes-table', join(dir, 'genes.tsv')]);
    const t = wasm.build_transcriptome(gunzipSync(readFileSync(FASTA_GZ)), gunzipSync(readFileSync(GFF_GZ)).toString());
    assert.ok(t.n_transcripts > 6000);
    assert.ok(Buffer.from(t.fasta).equals(readFileSync(join(dir, 't.fa'))), 'transcripts.fa differs');
    assert.equal(t.tx2gene, readFileSync(join(dir, 't2g.tsv'), 'utf8'), 'tx2gene differs');
    assert.equal(t.genes_tsv, readFileSync(join(dir, 'genes.tsv'), 'utf8'), 'genes table differs');
    assert.equal(t.gene_info_tsv, readFileSync(join(dir, 'gene_info.tsv'), 'utf8'), 'gene_info differs');
    t.free();
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

// Gene sums made without the Rust code: sum est_counts per gene, round once per gene.
function jsGeneCounts(abundance, tx2gene, scale) {
  const t2g = new Map(tx2gene.trim().split('\n').map((l) => l.split('\t')));
  const sums = new Map();
  const lines = abundance.trim().split('\n');
  const cols = lines[0].split('\t');
  const idc = cols.indexOf('target_id');
  const ec = cols.indexOf('est_counts');
  for (const l of lines.slice(1)) {
    const f = l.split('\t');
    const g = t2g.get(f[idc]);
    sums.set(g, (sums.get(g) || 0) + Number(f[ec]) * scale);
  }
  return new Map([...sums].map(([g, s]) => [g, Math.round(s)]));
}

test('gene_counts and assemble_matrix give the bytes of `qrnamap matrix`', { skip }, () => {
  const dir = mkdtempSync(join(tmpdir(), 'qrnamap-wasm-'));
  try {
    // The CLI makes tx2gene; both paths use it.
    run(CLI, ['transcripts', '--reference', FASTA_GZ, '--gff', GFF_GZ,
      '--out', join(dir, 't.fa'), '--tx2gene', join(dir, 't2g.tsv')]);
    const t2g = readFileSync(join(dir, 't2g.tsv'), 'utf8');
    const abundance = readFileSync(ABUNDANCE, 'utf8');
    // Three samples from the one native abundance table: as it is, again under a second name,
    // and with est_counts multiplied by 0.5 (so the sums differ and rounding is exercised).
    const scales = { S1: 1, S2: 1, S3: 0.5 };
    const scaled = (s) => {
      const lines = abundance.trim().split('\n');
      const ec = lines[0].split('\t').indexOf('est_counts');
      return [lines[0], ...lines.slice(1).map((l) => {
        const f = l.split('\t');
        f[ec] = String(Number(f[ec]) * s);
        return f.join('\t');
      })].join('\n') + '\n';
    };

    // WASM path.
    const samples = [];
    for (const [name, s] of Object.entries(scales)) {
      const gc = JSON.parse(wasm.gene_counts(scaled(s), t2g));
      samples.push({ name, counts: gc.counts });
      // Independent check of the sums.
      const js = jsGeneCounts(abundance, t2g, s);
      assert.deepEqual(Object.entries(gc.counts).sort(), [...js].sort(), `gene sums of ${name}`);
    }
    const m = JSON.parse(wasm.assemble_matrix(JSON.stringify({ tx2gene: t2g, samples })));

    // CLI path: per-sample gene count files, summary files, then `qrnamap matrix`.
    run('zstd', ['-q', '-o', join(dir, 't2g.tsv.zst'), join(dir, 't2g.tsv')]);
    let sheet = 'sample\tfastq1\n';
    for (const [name, s] of Object.entries(scales)) {
      const js = jsGeneCounts(abundance, t2g, s);
      const body = `gene\t${name}\n` + [...js].sort().map(([g, n]) => `${g}\t${n}\n`).join('');
      writeFileSync(join(dir, `${name}.genecounts.tsv`), body);
      run('zstd', ['-q', '-o', join(dir, `${name}.genecounts.tsv.zst`), join(dir, `${name}.genecounts.tsv`)]);
      const summary = 'sample\tmode\tquantifier\treads_processed\tfraction_mapped\tfraction_assigned\tstrandedness\tqrnamap_version\n' +
        `${name}\tquant\tkallisto\t1\t1\t1\tunstranded\t0.1.0\n`;
      writeFileSync(join(dir, `${name}.summary.tsv`), summary);
      run('zstd', ['-q', '-o', join(dir, `${name}.summary.tsv.zst`), join(dir, `${name}.summary.tsv`)]);
      sheet += `${name}\t${name}.fq\n`;
    }
    writeFileSync(join(dir, 'samples.tsv'), sheet);
    run(CLI, ['matrix', '--out', dir, '--samples', join(dir, 'samples.tsv'),
      '--genes-from', 'tx2gene', join(dir, 't2g.tsv.zst')]);
    const cli = run('zstd', ['-dc', join(dir, 'counts.tsv.zst')]).toString();
    assert.ok(cli.split('\n').length > 6000);
    assert.equal(m.tsv, cli, 'counts.tsv differs');
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test('qc, kallisto_args, de_command and validate_sheet', () => {
  const q = JSON.parse(wasm.qc('gene\tA\tB\tC\ng1\t10\t200\t5\ng2\t500\t20\t90\ng3\t40\t40\t40\n'));
  assert.equal(q.correlation.length, 3);
  assert.equal(q.pca.scores.length, 3);
  assert.throws(() => wasm.qc('gene\tA\tB\ng1\t1\t0\n'), /B/);
  assert.deepEqual(JSON.parse(wasm.kallisto_args('{"kind":"index","fasta":"t.fa","index":"t.idx"}')), ['index', '-i', 't.idx', 't.fa']);
  assert.match(wasm.de_command('s.tsv', 'c.tsv', '~ condition', 'condition,a,b'), /^qrnamap de --counts c.tsv --samples s.tsv /);
  assert.equal(JSON.parse(wasm.validate_sheet('sample\tfastq1\nA\ta.fq\nA\tb.fq\n')).ok, false);
  assert.throws(() => wasm.build_transcriptome(new Uint8Array(0), 'x'), /empty/);
});
