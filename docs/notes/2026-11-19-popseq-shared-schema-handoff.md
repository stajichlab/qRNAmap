# Note for PopSeq: a shared genome-assembly schema with `genome-catalog`

Date: 2026-11-19. From: qRNAmap (design in `docs/superpowers/specs/2026-11-19-qrnamap-design.md`). Status: proposal, no change to PopSeq is requested yet. DRAFT: not copied into the PopSeq repository.

## What this is

Three lab tools need a list of fungal genome assemblies with download addresses: qmutfind, qRNAmap and PopSeq. Today each builds its own copy from NCBI's `assembly_summary.txt`. We plan a small Rust crate, `genome-catalog`, that owns the schema and the derived layer (FASTA and GFF3 URLs, file sizes, annotation facts). It will start inside the qRNAmap workspace and later move to its own repository, `stajichlab/genome-catalog`.

We ask whether PopSeq's assembly table can be the upstream input, and whether the two schemas can share a key and column names.

## What exists now (read 2026-11-19)

PopSeq: `catalog/kingdom/fungi_assemblies.csv.gz`, built by `scripts/popseq_assemblies.py` (Python standard library only). 26,826 lines including the header (rows not counted separately). Columns:
`assembly_accession, refseq_accession, asm_name, sample_accession, study_accession, tax_id, species_taxid, species_name, organism_name, strain, isolate, assembly_level, genome_rep, refseq_category, relation_to_type_material, seq_rel_date, asm_submitter, genome_size, excluded_from_refseq`.
One row per GenBank assembly (GCA).

qmutfind: `site/genomes/genomes.json`, about 50 hand-chosen RefSeq genomes plus two GenBank ones. Per genome:
`accession, source, label, organism, strain, level, version_status, assembly, released, genome_size, genes, protein_coding_genes, annotation, fasta{name,url,bytes}, gff{name,url,bytes}`.
Sizes come from HEAD requests. A genome with no GFF3 has `gff: null`.

## Overlap and differences

| Concept | PopSeq column | qmutfind field |
|---|---|---|
| assembly id | `assembly_accession` (GCA), `refseq_accession` | `accession` (GCF or GCA) |
| assembly name | `asm_name` | `assembly` |
| organism, strain | `organism_name`, `strain` | `organism`, `strain` |
| level | `assembly_level` | `level` |
| release date | `seq_rel_date` | `released` |
| size | `genome_size` | `genome_size` |
| taxonomy, BioSample, study | `tax_id`, `species_taxid`, `sample_accession`, `study_accession`, `isolate` | none |
| download URLs and sizes in bytes | none that I saw | `fasta`, `gff` |
| gene counts, annotation source | none | `genes`, `protein_coding_genes`, `annotation` |

## Proposal

1. **Key:** `assembly_accession` (GCA) is the shared key. A RefSeq genome is linked by `refseq_accession`. qmutfind keys on the GCF accession for RefSeq entries today, so `genome-catalog` would store both.
2. **Direction:** PopSeq's table flows into `genome-catalog`. `genome-catalog` adds columns (URLs, byte sizes, annotation facts). It does not change PopSeq's columns or its build.
3. **Shared names:** where a concept is in both tables, use PopSeq's column name (`assembly_level`, `asm_name`, `genome_size`, `seq_rel_date`). qmutfind's JSON would be generated from the shared table, so its field names would change in a later qmutfind release.
4. **Format:** `genome-catalog` publishes a versioned file (gzip CSV or JSON, with a schema version field). Python tools read the file. They do not link the Rust code.
5. **Scope:** genome FASTA and GFF3 only. Reads, ENA run metadata and run selection stay in PopSeq.

## Questions for PopSeq's owner

- Can PopSeq's NCBI parsing in `popseq_assemblies.py` be the single parser, with `genome-catalog` calling its output? Or would you rather `genome-catalog` own the parser and PopSeq read the result? I lean to the first, because the script exists and is tested. I have not read its tests.
- Does the table contain only the current NCBI summary (suppressed assemblies are absent, per the script's header)? A fixed genome list for a published analysis may need suppressed ones. Should the shared file keep a dated snapshot?
- Are FASTA and GFF3 addresses derivable from `assembly_accession` and `asm_name`? The addresses in qmutfind's JSON follow the pattern `.../genomes/all/GCF/028/009/165/<accession>_<asm_name>/<accession>_<asm_name>_genomic.fna.gz`. I saw that pattern in one entry and did not check all of them, and I did not check that GFF3 files exist for every assembly (qmutfind notes one that has none).
- Where do gene counts and the annotation source come from for qmutfind's list? I have not traced that. If they come from a source PopSeq does not read, `genome-catalog` will fetch them.
- Does PopSeq's `catalog/index.json` need a new entry for a table it consumes but does not build?

## What is not decided

- The exact schema. It is frozen only after the columns above are compared on real rows.
- Ensembl. Its download behavior and CORS headers are untested, so the schema has no Ensembl fields yet.
- The release process for the shared file.

## Next step if you agree

A short comparison script in `genome-catalog` that reads PopSeq's CSV and qmutfind's JSON, joins them on accession, and reports rows that disagree. That gives us measured overlap before any schema is fixed.
