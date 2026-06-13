# AMR Resistance Gene Detection Module
### *Neisseria gonorrhoeae* Diagnostic Pipeline

> **Module:** Days 1–21 of the 60–90 day bioinformatics consulting portfolio  
> **Repository:** [bkhimek/N_gonorrhoeae-diagnostic-pipeline](https://github.com/bkhimek/N_gonorrhoeae-diagnostic-pipeline)  
> **Last updated:** June 2026

---

## Overview

This module adds antimicrobial resistance (AMR) profiling to the existing *N. gonorrhoeae* pan-genome and primer design pipeline. It processes all 283 assembled genomes through a two-pronged approach:

1. **AMRFinderPlus** — gene-level detection (β-lactamases, efflux genes, plasmid-mediated resistance)
2. **Custom alignment-based mutation calling** — clinical point mutations in penA, mtrR, gyrA, parC, and porB, anchored to the FA-1909 reference strain

---

## Target Genes and Resistance Markers

| Gene | Resistance mechanism | Key mutations tracked |
|------|---------------------|----------------------|
| **gyrA** | Fluoroquinolone (DNA gyrase) | S91F/Y, D95G/N/A/H |
| **parC** | Fluoroquinolone (topoisomerase IV) | D86N/Y, S87R/N/I/T, E91K/G |
| **penA** | β-Lactam (PBP2) | Mosaic allele (≥5 aa substitutions), A345G, D346ins, A501V, G543S/D |
| **mtrR** | Efflux pump (MtrCDE repressor) | A39T |
| **porB** | Porin permeability | G120D/K/N, A121D/P/N |

---

## Pipeline Architecture

```
prokka_output/{strain}/
     └── {strain}.faa  ──► BLAST vs FA-1909 reference genes
                                    │
                             per-strain target FASTAs
                                    │ (collect all 283 strains)
                             COLLECT_GENES
                                    │
                    ┌───────────────┼───────────────┐
                  gyrA            parC         penA / mtrR / porB
                    │               │               │
                  MAFFT           MAFFT           MAFFT
                    │               │               │
              parse_mutations  parse_mutations  parse_mutations
                    └───────────────┼───────────────┘
                                    │
               AMRFinderPlus ───────┤
               (protein-only)       │
                                    ▼
                           RESISTANCE_PROFILE
                                    │
                    ┌───────────────┴────────────────┐
            resistance_profiles.tsv          amr_summary.html
```

---

## Reference Strain

**FA-1909** (GCF_000006845.1, GenBank AE004969) is used as the alignment anchor for all mutation numbering. This is the canonical *N. gonorrhoeae* reference genome. All amino acid positions reported follow FA-1909 numbering.

> WHO reference strains (WHO F, K, etc.) can be added as validation controls but are **not** used as the alignment reference.

---

## Tools and Versions

| Tool | Version | Purpose |
|------|---------|---------|
| AMRFinderPlus | 4.0.23 | Gene-level AMR detection |
| AMRFinder DB | 2025-07-16.1 | NCBI AMR database |
| BLAST+ | system | Gene extraction from assemblies |
| MAFFT | ≥7.0 | Multiple sequence alignment |
| Nextflow | 26.04.3 | Pipeline orchestration (DSL2) |
| Python | ≥3.9 | Mutation calling and reporting |

---

## Directory Structure

```
n_gonorrhoeae_project/
├── main_amr.nf                    # Pipeline entry point
├── nextflow.config                # Parameters and resource config
├── modules/amr/
│   ├── amrfinderplus.nf
│   ├── extract_targets.nf
│   ├── collect_genes.nf
│   ├── align_targets.nf
│   ├── call_mutations.nf
│   └── resistance_profile.nf
├── subworkflows/
│   └── amr_analysis.nf
├── bin/
│   ├── parse_mutations.py         # Alignment-based mutation caller
│   └── summarize_amr.py          # Profile merger and HTML reporter
├── reference/
│   └── genes/
│       ├── gyrA.faa               # FA-1909 reference sequences
│       ├── parC.faa
│       ├── penA.faa
│       ├── mtrR.faa
│       └── porB.faa
└── results/amr/
    ├── amrfinderplus/             # Per-strain AMRFinder TSVs
    ├── blast_targets/             # Per-strain extracted gene FASTAs
    ├── alignments/                # Per-gene MAFFT alignments
    ├── mutations/                 # Per-gene mutation tables
    ├── profiles/                  # resistance_profiles.tsv
    └── reports/                   # amr_summary.html
```

---

## Setup

### 1. Install dependencies
```bash
conda install -y -c bioconda ncbi-amrfinderplus mafft blast samtools
```

### 2. Update AMRFinder database
```bash
amrfinder -u
# Database installs to: ~/miniconda3/share/amrfinderplus/data/latest
```

### 3. Extract FA-1909 reference gene sequences (run once)
```bash
FAA=prokka_output/GCF_000006845.1/GCF_000006845.1.faa
samtools faidx $FAA

samtools faidx $FAA BMGINFML_01541 > reference/genes/penA.faa
samtools faidx $FAA BMGINFML_01361 > reference/genes/mtrR.faa
samtools faidx $FAA BMGINFML_00633 > reference/genes/gyrA.faa
samtools faidx $FAA BMGINFML_01250 > reference/genes/parC.faa
samtools faidx $FAA BMGINFML_01819 > reference/genes/porB.faa

for gene in penA mtrR gyrA parC porB; do
    sed -i "1s/.*/>${gene}_FA1909/" reference/genes/${gene}.faa
done
```

---

## Running the Pipeline

```bash
cd /home/himek/n_gonorrhoeae_project

# Full run — all 283 strains (~30–40 min on 8 cores)
nextflow run main_amr.nf -profile local

# Resume after interruption
nextflow run main_amr.nf -profile local -resume
```

---

## Output Files

### `results/amr/profiles/resistance_profiles.tsv`
One row per strain with 17 tab-separated columns:

| Column | Description |
|--------|-------------|
| strain_id | GCF accession |
| GyrA_S91 | S91F/Y mutation (True/False) |
| GyrA_D95 | D95G/N/A/H mutation |
| ParC_D86 | D86N/Y mutation |
| ParC_S87 | S87R/N/I/T mutation |
| ParC_E91 | E91K/G mutation |
| PenA_mosaic | ≥5 aa substitutions vs FA-1909 |
| PenA_A345G … PenA_G543 | Individual penA mutations |
| MtrR_A39T | Efflux regulator mutation |
| PorB_G120, PorB_A121 | Porin mutations |
| amrfinder_hits | Pipe-separated AMRFinder gene hits |
| predicted_cipro_r | GyrA S91 or D95 mutation present |
| predicted_ceph_r | penA mosaic allele present |

### `results/amr/reports/amr_summary.html`
Colour-coded summary table: green = wild-type, red = resistance mutation detected.

---

## Key Findings — 283 Strains (June 2026 Run)

| Resistance marker | N | Prevalence |
|-------------------|---|------------|
| GyrA S91F/Y | 169 | 60% |
| GyrA D95G/N/A | 158 | 56% |
| ParC S87R/N/I | 99 | 35% |
| penA mosaic allele | 143 | 51% |
| MtrR A39T | 71 | 25% |
| **Predicted ciprofloxacin R** | **170** | **60%** |
| **Predicted reduced cephalosporin susceptibility** | **143** | **51%** |

Notable AMRFinder-only findings: `blaTEM-1` (PPNG — plasmid-mediated penicillinase) and `tet(M)` (tetracycline resistance) detected in select strains — genes not detectable by point-mutation calling alone, demonstrating the value of the dual-approach design.

---

## Design Decisions

**Why protein-only AMRFinderPlus mode?**  
Prokka GFF3 uses `ID=`/`locus_tag=` attributes; AMRFinderPlus expects `protein_id=` (GenBank style). These are incompatible. Protein-only mode detects resistance genes and protein-level point mutations (POINTP) with sufficient sensitivity for the cross-check role in this pipeline.

**Why custom mutation calling instead of relying on AMRFinderPlus?**  
AMRFinderPlus with `--organism Neisseria_gonorrhoeae` is curated primarily for *N. meningitidis*. For QRDR point mutations specifically relevant to gonorrhoea treatment (GyrA S91, D95; ParC S87), alignment to FA-1909 with documented gonococcal mutation catalogues provides more reliable and interpretable results.

---

## Citation

- Feldgarden et al. (2021) AMRFinderPlus. *Scientific Reports* 11:12115  
- Katoh & Standley (2013) MAFFT. *Mol Biol Evol* 30:772–780  
- Seemann (2014) Prokka. *Bioinformatics* 30:2068–2069
