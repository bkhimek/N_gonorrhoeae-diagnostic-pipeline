# N. gonorrhoeae Diagnostic Target Discovery Pipeline

A Nextflow DSL2 pipeline for discovering species-specific PCR diagnostic targets in bacterial pathogens, applied here to *Neisseria gonorrhoeae*.

## Overview

The pipeline performs:
1. **Gene annotation** — Prokka on all input genomes
2. **Pan-genome analysis** — Panaroo identifies core genes (≥95% prevalence)
3. **Exclusivity screening** — BLAST core genes vs. nt_prok, excluding target species taxids
4. **Coverage analysis** — Union-of-hits interval merging to identify genuinely unique regions
5. **Primer design** — Primer3 on extracted unique region sequences

## Requirements

- Nextflow ≥ 22.10
- Docker (pipeline runs all tools in containers)
- NCBI nt_prok BLAST database (local copy)
- ~50 GB disk space for a run of ~300 genomes

## Quick Start

```bash
# Clone repository
git clone https://github.com/<your-username>/ngono-diagnostic-pipeline.git
cd ngono-diagnostic-pipeline

# Run with N. gonorrhoeae config
nextflow run main.nf \
    -profile docker \
    -c conf/n_gonorrhoeae.config \
    --input /path/to/genomes/ \
    --outdir results/
```

## Directory Structure

```
.
├── main.nf                        # Main pipeline
├── nextflow.config                # Docker, resource, volume config
├── conf/
│   └── n_gonorrhoeae.config       # Species-specific parameters
├── modules/
│   ├── prokka.nf
│   ├── panaroo.nf
│   ├── blast_exclusivity.nf
│   └── ...
├── results/
│   ├── panaroo/
│   │   └── combined_DNA_CDS.fasta # Pan-genome nucleotide sequences
│   ├── exclusivity/
│   │   ├── blast_results.tsv      # Full BLAST output (18 columns)
│   │   └── exclusivity.tsv        # PASS/FAIL per gene
│   ├── pass_fastas/
│   │   └── pass/                  # Per-gene FASTA for PASS genes
│   └── inclusivity/
├── scripts/
│   └── unique_regions.py          # Post-pipeline coverage analysis
└── unique_regions_for_primers.fasta  # Final primer design input
```

## Configuration

### `conf/n_gonorrhoeae.config`

Key parameters:

```groovy
params {
    species              = "Neisseria_gonorrhoeae"
    taxid                = 485
    negative_taxidlist   = "/path/to/gonorrhoeae_all_taxids.txt"  // 90 taxids
    blast_db             = "/path/to/nt_prok"
    blast_evalue         = "1e-5"
    core_threshold       = 0.95   // ≥95% genomes = core gene
    min_unique_region    = 150    // bp
}
```

### Negative taxid list

The `negative_taxidlist` file must contain all *N. gonorrhoeae* taxids (species + all strain-level children) to exclude them from BLAST. Generate with:

```bash
# Get all N. gonorrhoeae strain taxids from NCBI Taxonomy
# taxid 485 = N. gonorrhoeae species node
efetch -db taxonomy -id 485 -format xml | \
    xtract -pattern Taxon -element TaxId > gonorrhoeae_all_taxids.txt
```

The pipeline uses BLAST `-negative_taxidlist` rather than `-negative_taxids` to support large exclusion lists (>100 taxids).

## BLAST Output Columns

`blast_results.tsv` contains 18 tab-separated columns:

| # | Field | Description |
|---|-------|-------------|
| 1 | qseqid | Query gene name |
| 2 | sseqid | Subject accession |
| 3 | pident | % identity |
| 4 | length | Alignment length |
| 5 | mismatch | Mismatches |
| 6 | gapopen | Gap openings |
| 7 | qstart | Query start |
| 8 | qend | Query end |
| 9 | sstart | Subject start |
| 10 | send | Subject end |
| 11 | evalue | E-value |
| 12 | bitscore | Bit score |
| 13 | qlen | Query length |
| 14 | slen | Subject length |
| 15 | qcovs | Query coverage (%) |
| 16 | sscinames | Subject scientific name |
| 17 | staxids | Subject taxids |
| 18 | stitle | Subject title |

## Post-Pipeline Analysis: Unique Region Discovery

After the pipeline, run `scripts/unique_regions.py` to identify genuinely uncovered gene regions via interval merging:

```bash
python3 scripts/unique_regions.py
```

This merges all BLAST hit intervals per gene and finds gaps — positions not covered by any non-gonococcal organism. Output:

```
Gene       GeneLen  5prime_uncov  3prime_uncov  internal_gaps  max_unique_bp  pos
polA          2823           910             0              0            910  1-910
norR          1727             0           713              0            713  1015-1727
pyrI           525           161             0              0            161  1-161
bcp            435             1             0              0              1  1-1
...
```

## Results Summary (N. gonorrhoeae)

Input: 283 complete RefSeq genomes | Excluded: 90 *N. gonorrhoeae* taxids

| Gene | Unique Region | Size | Recommendation |
|------|---------------|------|----------------|
| polA | pos 1–910 (5′) | 910 bp | **Primary PCR target** |
| norR | pos 1015–1727 (3′) | 713 bp | **Primary PCR target** |
| pyrI | pos 1–161 (5′) | 161 bp | Marginal — needs validation |

See `case_study_neisseria_gonorrhoeae.md` for full analysis and biological interpretation.

## Primer Design

Primers were designed with Primer3 using sequences from `unique_regions_for_primers.fasta`:

```bash
# Primer3 parameters
PRIMER_OPT_SIZE=20 / MIN=18 / MAX=25
PRIMER_OPT_TM=60.0 / MIN=57.0 / MAX=63.0
PRIMER_MIN_GC=40.0 / MAX=65.0
PRIMER_PRODUCT_SIZE_RANGE=150-400
```

Candidates should be validated in-silico with NCBI Primer-BLAST (taxid 485) before wet-lab testing.

## Citation

If using this pipeline, please cite the underlying tools:
- **Prokka**: Seemann T. (2014) Bioinformatics 30:2068-2069
- **Panaroo**: Tonkin-Hill G. et al. (2020) Genome Biology 21:180
- **BLAST+**: Camacho C. et al. (2009) BMC Bioinformatics 10:421
- **Primer3**: Untergasser A. et al. (2012) Nucleic Acids Res 40:e115
