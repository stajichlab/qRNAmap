# Test sets for the browser page

The page offers two public test sets (Task 8 of `docs/superpowers/plans/2026-11-22-browser-page-v1.md`), defined in `site/testsets.json`:

1. **Yeast amphotericin B (GSE80357)**:
   - Organism: *Saccharomyces cerevisiae* (RefSeq GCF_000146045.2)
   - GEO series: GSE80357 / BioProject PRJNA318684
   - 3+3 matched paired-end HiSeq 2000 runs: Control (SRR3396385, SRR3396386, SRR3396387) vs AmpB (SRR3396381, SRR3396382, SRR3396384)
   - Reads: 101 bp paired-end, cut at N = 200,000 pairs
   - Contrast: `condition,AmpB,control`

2. ***Cryptococcus neoformans* H99 Lkh1 (PRJNA976751)**:
   - Organism: *Cryptococcus neoformans* var. *grubii* H99 (RefSeq GCF_000149245.1)
   - BioProject: PRJNA976751 / GEO GSE233612
   - Paper: Front Cell Infect Microbiol 2024;14:1367098 (PMID 38774630, PMC11106425)
   - 3+3 basal condition runs: WT (SRR24750781, SRR24750780, SRR24750777) vs *lkh1* deletion (SRR24750779, SRR24750774, SRR24750773)
   - Reads: 150 bp paired-end, cut at N = 75,000 pairs
   - Contrast: `genotype,lkh1_delta,wild_type`

## Fetch method

The browser streams the FASTQ files over HTTPS from ENA (`https://ftp.sra.ebi.ac.uk/vol1/fastq/...`), decompresses them as a stream with `gz.js` (fflate `VerifyingGunzip`), cuts at N records per mate, checks that record counts match and read names match, and gives plain FASTQ to kallisto in the worker.
In tests, the network requests are served by the local static server to avoid reliance on external ENA availability.
