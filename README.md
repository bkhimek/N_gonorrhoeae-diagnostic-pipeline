# Diagnostic Inclusivity / Exclusivity Pipeline

A reproducible, containerised Nextflow pipeline for in-silico validation of molecular diagnostic targets — ready for FDA/CE-IVDR regulatory submissions.

## What it does

Given a set of annotated genomes, this pipeline:

1. Builds a pan-genome (Panaroo) across all target species strains
2. Identifies conserved core genes present in ≥99% of strains with ≥98% sequence identity (inclusivity)
3. BLASTs candidate sequences against a panel of near-neighbour species (exclusivity)
4. Outputs a ranked list of diagnostic targets with PASS/REJECT decisions and nearest phylogenetic neighbour identity

## Why it matters

Regulatory submissions for molecular diagnostics (qPCR, isothermal assays) require documented evidence that the target sequence is:
- **Inclusive**: present and conserved across the target species
- **Exclusive**: absent or sufficiently divergent in related species

This pipeline automates that evidence generation in a fully traceable, reproducible workflow.

## Quick start

```bash
# Clone the repo
git clone https://github.com/bkhimek/diagnostic-inclusivity-exclusivity-pipeline.git
cd diagnostic-inclusivity-exclusivity-pipeline

# Run with your own data
nextflow run main.nf \
  --gffs "/path/to/prokka_output/*/*.gff" \
  --blast_db "/path/to/non_target_species_db" \
  --core_threshold 0.99 \
  --identity_thr 98.0 \
  --min_pident 85.0 \
  --min_qcovs 80.0
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--gffs` | required | Glob path to Prokka GFF files |
| `--blast_db` | required | Path to BLAST database of non-target species |
| `--core_threshold` | 0.99 | Minimum fraction of strains gene must be present in |
| `--identity_thr` | 98.0 | Minimum average pairwise identity (%) for inclusivity |
| `--min_pident` | 85.0 | Identity threshold for exclusivity REJECT decision |
| `--min_qcovs` | 80.0 | Coverage threshold for exclusivity REJECT decision |

## Requirements

- Nextflow ≥ 26.0
- Docker

All tools (Panaroo, BLAST, MAFFT, Python scripts) run in containers — no conda environments needed.

## Output

| File | Description |
|---|---|
| `results/inclusivity/candidates.tsv` | Genes passing inclusivity filter |
| `results/exclusivity/exclusivity.tsv` | Per-gene exclusivity decisions with nearest neighbour |
| `results/pass_fastas/pass/*.fasta` | Final PASS candidate sequences |

## Case study: Streptococcus pyogenes (Strep A)

Applied to 375 complete S. pyogenes genomes against a panel of 15 near-neighbour streptococcal species:
- 494 core genes identified at ≥99% presence, ≥98% identity
- 463 genes passed exclusivity screening
- Top candidates show >25% sequence divergence from nearest neighbour

## Citation

If you use this pipeline, please cite:
- [Panaroo](https://github.com/gtonkinhill/panaroo)
- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK153387/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
