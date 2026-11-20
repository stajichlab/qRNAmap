use genome_catalog::{builtin, parse, size_warning, CatalogError, FileRef, Genome};

fn genome_of_size(size: u64) -> Genome {
    Genome {
        accession: "GCF_TEST.1".into(),
        source: "RefSeq".into(),
        label: "Testus fungus".into(),
        organism: "Testus fungus".into(),
        strain: None,
        level: "Chromosome".into(),
        assembly: "ASM1".into(),
        genome_size: size,
        genes: None,
        fasta: FileRef { name: "a.fna.gz".into(), url: "https://x/a.fna.gz".into(), bytes: 1 },
        gff: None,
    }
}

#[test]
fn builtin_parses() {
    let c = builtin();
    assert_eq!(c.schema_version, 1);
    assert!(!c.generated.is_empty());
    assert!(c.genomes.len() > 100);
}

#[test]
fn find_known_accession() {
    let c = builtin();
    let g = c.find("GCF_000002655.1").expect("Af293 is in the catalog");
    assert_eq!(g.organism, "Aspergillus fumigatus Af293");
    assert_eq!(g.strain.as_deref(), Some("Af293"));
    assert_eq!(g.genome_size, 29_384_958);
    assert_eq!(g.genes, Some(10086));
    assert_eq!(
        g.fasta.url,
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz"
    );
    assert_eq!(g.fasta.bytes, 9_003_172);
    assert!(c.find("GCF_NOPE.1").is_none());
}

#[test]
fn genome_without_gff_not_usable_for_rna() {
    let c = builtin();
    let cea10 = c.find("GCA_051225625.1").expect("CEA10 is in the catalog");
    assert!(cea10.gff.is_none());
    let usable = c.usable_for_rna();
    assert!(usable.iter().all(|g| g.gff.is_some()));
    assert!(!usable.iter().any(|g| g.accession == "GCA_051225625.1"));
    assert!(usable.iter().any(|g| g.accession == "GCF_000002655.1"));
    assert_eq!(usable.len(), c.genomes.len() - c.genomes.iter().filter(|g| g.gff.is_none()).count());
}

#[test]
fn size_warning_above_50mb() {
    assert_eq!(size_warning(&genome_of_size(50_000_000)), None);
    assert_eq!(size_warning(&genome_of_size(29_384_958)), None);
    let w = size_warning(&genome_of_size(67_000_000)).expect("67 Mb is above the limit");
    assert!(w.contains("Testus fungus"), "{w}");
    assert!(w.contains("67.0 Mb"), "{w}");
    assert!(size_warning(&genome_of_size(50_000_001)).is_some());
}

#[test]
fn unknown_schema_version_is_error() {
    let json = r#"{"schema_version": 2, "generated": "2026-01-01", "genomes": []}"#;
    assert_eq!(parse(json).unwrap_err(), CatalogError::UnknownSchemaVersion(2));
}

#[test]
fn missing_schema_version_is_error() {
    let json = r#"{"generated": "2026-01-01", "genomes": []}"#;
    assert!(matches!(parse(json), Err(CatalogError::Json(_))));
}

#[test]
fn malformed_json_is_error() {
    assert!(matches!(parse("{not json"), Err(CatalogError::Json(_))));
}

#[test]
fn null_or_absent_optional_fields_and_unknown_fields_are_accepted() {
    let json = r#"{"schema_version": 1, "generated": "g", "extra": 1, "genomes": [
        {"accession": "A", "source": "s", "label": "l", "organism": "o", "level": "L", "assembly": "a",
         "genome_size": 5, "strain": null, "genes": null, "gff": null, "annotation": "x",
         "fasta": {"name": "n", "url": "u", "bytes": 3}}]}"#;
    let c = parse(json).unwrap();
    assert_eq!(c.genomes[0].strain, None);
    assert_eq!(c.genomes[0].genes, None);
    assert!(c.genomes[0].gff.is_none());
}
