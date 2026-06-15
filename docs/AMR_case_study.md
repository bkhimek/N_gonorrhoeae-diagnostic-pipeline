# Bioinformatics Case Study: Antimicrobial Resistance Profiling at Scale

**A Dual-Approach AMR Module for 283 *Neisseria gonorrhoeae* Genomes**

> Pipeline repository: [bkhimek/N_gonorrhoeae-diagnostic-pipeline](https://github.com/bkhimek/N_gonorrhoeae-diagnostic-pipeline)
> June 2026

---

## Executive Summary

Gonorrhoea, caused by *Neisseria gonorrhoeae*, is a major public health concern because *N. gonorrhoeae* has acquired resistance to nearly every antibiotic class used for treatment, making genomic surveillance a public health priority.

The module processed 283 assembled genomes through two complementary approaches: AMRFinderPlus for gene-level detection of acquired resistance determinants, and a custom alignment-based pipeline for clinical point mutations in five key genes (*penA*, *gyrA*, *parC*, *mtrR*, *porB*) anchored to the FA-1909 reference strain.

| Parameter | Value |
|---|---|
| Genomes analysed | 283 *N. gonorrhoeae* assemblies |
| Reference strain | FA-1909 (GCF_000006845.1) |
| Pipeline framework | Nextflow DSL2, version 26.04.3 |
| Gene-level tool | AMRFinderPlus 4.0.23, database 2025-07-16.1 |
| Predicted cipro resistance | 170/283 strains (60%) |
| Predicted reduced ceph susceptibility | 143/283 strains (51%) |

---

## 1. Background and Clinical Context

### 1.1 The AMR challenge in gonorrhoea

*N. gonorrhoeae* has a remarkable capacity to acquire and maintain antibiotic resistance. Penicillin resistance emerged in the 1970s (plasmid-encoded beta-lactamase), tetracycline resistance followed in the 1980s, and fluoroquinolone resistance spread globally through the 1990s-2000s via chromosomal mutations in DNA gyrase (GyrA) and topoisomerase IV (ParC). Third-generation cephalosporins became the last reliable monotherapy option, but strains with reduced susceptibility driven by mosaic *penA* alleles are now well documented.

### 1.2 Resistance mechanisms targeted

| Gene | Drug class | Mechanism |
|---|---|---|
| *gyrA* | Fluoroquinolone | QRDR mutations (S91, D95) in DNA gyrase subunit A |
| *parC* | Fluoroquinolone | QRDR mutations (D86, S87, E91) in topoisomerase IV |
| *penA* | Beta-lactam | Mosaic PBP2 (>=5 aa substitutions vs FA-1909) |
| *mtrR* | Multiple | A39T derepresses MtrCDE efflux pump |
| *porB* | Beta-lactam | G120/A121 reduce outer membrane permeability |

### 1.3 Reference strain rationale

FA-1909 (GCF_000006845.1) is the canonical, fully closed *N. gonorrhoeae* reference genome used in the majority of published AMR studies. Using FA-1909 numbering ensures positions reported here correspond directly to positions cited in clinical literature.

---

## 2. Methods

### 2.1 Pipeline architecture

**Track 1 - Gene-level detection (AMRFinderPlus)**
- `AMRFINDERPLUS`: each strain's protein FASTA passed to `amrfinder` in protein-only mode with `--organism Neisseria_gonorrhoeae`

**Track 2 - Point mutation calling (custom alignment pipeline)**
- `EXTRACT_TARGETS`: BLAST searches each proteome against FA-1909 reference gene sequences; hits passing 80% identity / 80% query coverage extracted with `samtools faidx`
- `COLLECT_GENES`: sequences from all 283 strains pooled per gene
- `ALIGN_TARGETS`: MAFFT multiple sequence alignment per gene
- `CALL_MUTATIONS`: `parse_mutations.py` compares aligned sequences to FA-1909 at all catalogued resistance positions; QRDR window scanned for novel variants
- `RESISTANCE_PROFILE`: `summarize_amr.py` merges outputs into `resistance_profiles.tsv` and `amr_summary.html`

### 2.2 Key design decisions

**Protein-only AMRFinderPlus mode**
Prokka GFF3 uses `ID=`/`locus_tag=` attributes; AMRFinderPlus expects `protein_id=` (GenBank convention). These formats are incompatible. Protein-only mode provides equivalent gene-level detection without the format dependency.

**Custom alignment over AMRFinderPlus mutation rules**
AMRFinderPlus `--organism Neisseria_gonorrhoeae` curates point mutation rules primarily for *N. meningitidis*. Alignment to FA-1909 with literature-validated mutation catalogues provides more reliable results for gonorrhoea-specific QRDR positions. The two-track approach is complementary: AMRFinderPlus captures acquired genes (blaTEM-1, tet(M)) invisible to the alignment pipeline.

### 2.3 Software versions

| Tool | Version | Purpose |
|---|---|---|
| Nextflow | 26.04.3 | Pipeline orchestration (DSL2) |
| AMRFinderPlus | 4.0.23 | Gene-level AMR detection |
| AMRFinder database | 2025-07-16.1 | NCBI curated AMR gene database |
| BLAST+ | system | Gene extraction via protein search |
| MAFFT | >=7.0 | Multiple sequence protein alignment |
| Python | >=3.9 | Mutation calling and reporting |

---

## 3. Results

### 3.1 Resistance mutation prevalence (custom alignment pipeline)

| Marker | Mechanism | N strains | Prevalence |
|---|---|---|---|
| GyrA S91F/Y | Fluoroquinolone (DNA gyrase) | 169 | 60% |
| GyrA D95G/N/A | Fluoroquinolone (DNA gyrase) | 158 | 56% |
| ParC S87R/N/I/T | Fluoroquinolone (topoisomerase IV) | 99 | 35% |
| PenA mosaic allele | Cephalosporin (altered PBP2) | 143 | 51% |
| MtrR A39T | Efflux pump (MtrCDE) | 71 | 25% |

### 3.2 Predicted clinical resistance

| Predicted phenotype | Definition | N strains | Prevalence |
|---|---|---|---|
| Ciprofloxacin resistant | GyrA S91 OR GyrA D95 mutation | 170 | 60% |
| Reduced cephalosporin susceptibility | *penA* mosaic allele present | 143 | 51% |

### 3.3 AMRFinderPlus gene-level findings

AMRFinderPlus detected resistance-associated markers spanning six antibiotic classes, including determinants invisible to the alignment pipeline.

| Element | Class | Mechanism | N | % |
|---|---|---|---|---|
| gyrA_S91F | Fluoroquinolone | DNA gyrase QRDR | 169 | 60% |
| gyrA_D95A | Fluoroquinolone | DNA gyrase QRDR | 81 | 29% |
| gyrA_D95G | Fluoroquinolone | DNA gyrase QRDR | 62 | 22% |
| gyrA_D95N | Fluoroquinolone | DNA gyrase QRDR | 15 | 5% |
| parC_S87R | Fluoroquinolone | Topoisomerase IV QRDR | 80 | 28% |
| parC_D86N | Fluoroquinolone | Topoisomerase IV QRDR | 32 | 11% |
| parC_S87N | Fluoroquinolone | Topoisomerase IV QRDR | 16 | 6% |
| parC_E91G | Fluoroquinolone | Topoisomerase IV QRDR | 12 | 4% |
| penA_A510V | Beta-lactam | PBP2 mosaic substitution | 222 | 78% |
| penA_F504L | Beta-lactam | PBP2 mosaic substitution | 222 | 78% |
| penA_D346DD | Beta-lactam | PBP2 insertion | 176 | 62% |
| penA_A516G | Beta-lactam | PBP2 mosaic substitution | 141 | 50% |
| ponA_L421P | Beta-lactam | PBP1 mutation (synergistic) | 160 | 57% |
| pbp2 | Beta-lactam | Mosaic PBP2 allele (HMM) | 69 | 24% |
| blaTEM-1 | Beta-lactam | Plasmid beta-lactamase (PPNG) | 42 | 15% |
| blaTEM-135 | Beta-lactam | Extended-spectrum BL variant | 17 | 6% |
| mtrR_A39T | Multiple | MtrCDE efflux derepression | 77 | 27% |
| mtrR_G45D | Multiple | MtrCDE efflux derepression | 40 | 14% |
| porB1b_G120K | Beta-lactam | Porin permeability reduction | 107 | 38% |
| porB1b_A121D | Beta-lactam | Porin permeability reduction | 85 | 30% |
| rpsJ_V57M | Tetracycline | 30S ribosomal S10 mutation (chr.) | 211 | 75% |
| tet(M) | Tetracycline | Plasmid ribosomal protection (TRNG) | 48 | 17% |
| folP_R228S | Sulfonamide | DHPS active-site mutation | 201 | 71% |
| rpoB_H553N | Rifampicin | RNA polymerase beta subunit | 86 | 30% |

**Key observations:**

- **rpsJ_V57M (75%) vs tet(M) (17%):** chromosomal tetracycline resistance far outpaces plasmid-mediated TRNG, reflecting early selection pressure predating tet(M) emergence.
- **folP_R228S (71%):** sulfonamide resistance near-ubiquitous despite >50 years since sulfonamides were abandoned — classic resistance fixation.
- **ponA_L421P (57%):** PBP1 mutation acts synergistically with mosaic *penA*; cephalosporin resistance is polygenic in *N. gonorrhoeae*.
- **rpoB_H553N (30%):** rifampicin resistance common despite no clinical use, likely reflecting co-selection or clonal expansion.
- **blaTEM-135 (6%):** extended-spectrum BL variant warranting ongoing surveillance.
- **Dual-method value:** tet(M), blaTEM-1, rpsJ, folP, rpoB entirely invisible to alignment pipeline — demonstrating why the two-track design is complementary rather than redundant.

---

## 3.6 Core-genome phylogenetic analysis

Core-genome phylogenetics was performed using Parsnp v2.1.5, with FA1090 (GCF_000006845.1) as the reference. All 284 assemblies (283 clinical strains + reference) were included. Parsnp identified the shared core genome and generated a maximum-likelihood SNP tree; branch lengths represent SNP distance from the root. An interactive visualisation with AMR burden overlay is available in the dashboard Phylogeny tab.

### 3.6.1 Phylogenetic structure

The core-genome tree resolves the 283 clinical strains into **four major clades**. Three clades are compact and internally consistent, suggesting recent clonal expansion or founder effects. A fourth clade is more divergent, contributing greater within-group SNP diversity. The overall topology indicates that the 283-strain collection captures multiple independent lineages circulating within the sampled population, rather than a single outbreak clone.

### 3.6.2 AMR burden across the phylogeny

AMR burden — the count of unique resistance elements per strain as detected by AMRFinderPlus — was overlaid on the phylogeny as a colour scale (white = 0 elements, red = maximum):

| Metric | Value |
|--------|-------|
| Strains with AMR data | 283 |
| Minimum AMR burden | 0 elements |
| Maximum AMR burden | 23 elements |
| Mean AMR burden | 10.6 elements |
| Median AMR burden | 12 elements |

High AMR burden is not uniformly distributed across the tree. At least one clade is enriched for high-burden strains (≥ 15 resistance elements), consistent with clonal spread of a multidrug-resistant lineage. However, elevated AMR burden is also observed across multiple clades, indicating that some resistance determinants — particularly *gyrA* QRDR mutations (60% prevalence) and *rpsJ_V57M* (75% prevalence) — have been acquired convergently in independent lineages.

### 3.6.3 Clonality versus convergence

The co-existence of clade-specific resistance clusters and cross-clade high-prevalence markers has the following implications:

- **Clonal spread:** Resistance lineages in the compact clades likely descend from a common ancestor that had already acquired a core set of resistance mutations. Surveillance should monitor these lineages for additional resistance acquisition.
- **Convergent evolution:** The ubiquity of sulfonamide resistance (*folP_R228S*, 71%) and tetracycline resistance (*rpsJ_V57M*, 75%) across all clades confirms that these markers have become effectively fixed in the population through repeated independent selection, irrespective of phylogenetic background.
- **Actionable finding:** Strains carrying both fluoroquinolone (*gyrA*/*parC*) and cephalosporin (*penA* mosaic) resistance markers are a priority concern. Their distribution across clades suggests that treatment failure risk is not confined to a single emergent lineage.

---

## 3.7 Promoter-region variant analysis

In addition to coding-sequence mutations, *N. gonorrhoeae* resistance is modulated by regulatory mutations in gene promoters. The most clinically relevant is a 13-bp inverted repeat (IR; 5′-TAATTTTTATTTT-3′) in the *mtrR/mtrCDE* intergenic region: deletion of this IR derepresses the MtrCDE efflux pump, conferring reduced susceptibility to penicillins, tetracycline, and macrolides (Hagman & Shafer, 1995; PMID 7506352). A custom BLAST + MAFFT pipeline extracted and compared the 400 bp upstream of *mtrR* and 300 bp upstream of *porB* in all 283 assemblies against the FA1090 reference regions.

### 3.7.1 Reference genome caveat

FA1090 itself lacks the classic 13-bp IR in its *mtrR* upstream region. Consequently, clinical strains that **retain** the intact IR appear as insertions relative to FA1090, rather than the expected deletions. All indel counts are expressed as differences from the FA1090 promoter state, not from a susceptible-strain baseline.

### 3.7.2 mtrR upstream region

All 283 assemblies produced a BLAST hit (100% detection rate). Indel counts relative to FA1090:

| Indels vs FA1090 | Strains | % cohort | Interpretation |
|:---:|---:|---:|---|
| 0 | 1 | 0.4 | FA1090 itself (identical) |
| 1 | 157 | 55.5 | Minor structural variant (most common) |
| 2 | 112 | 39.6 | Two-position change |
| 3–9 | 8 | 2.8 | Moderate variant |
| 27 | 1 | 0.4 | ⚠ Possible assembly issue |
| 189 | 4 | 1.4 | ⚠ Likely assembly artefact |

The dominant pattern (1–2 indels, 269/283 = 95.1%) reflects either normal sequence polymorphism in the AT-rich promoter or IR retention in the majority of strains. Five strains with >20 indels were flagged; cross-referencing confirmed their AMR burden (8 elements each) is within the normal cohort range (0–23), indicating the issue is localised to this locus rather than reflecting genome-wide quality problems.

### 3.7.3 porB upstream region

The *porB* upstream region is highly conserved: 276/283 strains (97.5%) show no indels relative to FA1090; four carry a single indel; three show 160 indels (the same assembly-quality outliers above). No specific promoter variants were called in this region.

### 3.7.4 Limitations and next steps

- **Reference strain:** A susceptible FA19-type reference (with the IR intact) would enable direct IR-deletion calling and distinguish resistant from susceptible promoter types.
- **Validation:** Short-read mapping to the flagged strains would confirm whether the high indel counts reflect true structural variants or assembly errors.
- **Integration:** Combining *mtrR* promoter status with *mtrR_A39T* coding mutation data (27% prevalence in this cohort) would provide a more complete picture of MtrCDE efflux regulation across the 283 strains.

---
## 4. Technical Challenges and Solutions

| Challenge | Root cause | Resolution |
|---|---|---|
| AMRFinderPlus GFF mismatch | Prokka GFF3 incompatible with AMRFinderPlus format | Protein-only mode |
| BLAST silent no-hits | Dots in DB name + `-quiet` hid errors | Safe DB name `strain_db`; remove `-quiet` during dev |
| samtools faidx no output | Missing `.fai` index | Add `samtools faidx $faa` before extraction |
| `fromPath` empty channel | Nextflow 26 requires `type:'dir'` for directory globs | Add `type: 'dir'` |
| Processes show `-` | `projectDir` resolved to `/tmp/` | Always run from project root |
| `RESISTANCE_PROFILE` arg error | Collected files passed as space-separated list | `mkdir` + `mv` into subdir |

See [AMR_challenges_and_tips.md](AMR_challenges_and_tips.md) for the full troubleshooting reference.

---

## 5. Portfolio Narrative

### Skills demonstrated

- **Nextflow DSL2:** modules, subworkflows, channel operations, resume capability across 283 samples
- **Tool integration:** AMRFinderPlus, BLAST+, MAFFT, samtools with format-specific configuration
- **Python scripting:** custom mutation caller and reporting script
- **Problem diagnosis:** systematic debugging of five Nextflow-specific and three tool-specific issues
- **Biological interpretation:** results contextualised against published surveillance literature

### Potential extensions

- WHO reference strain panel (WHO F-Q) as validation controls
- Promoter-region variants via nucleotide alignment track
- Core-genome SNP phylogeny to map resistance onto strain clades
- Interactive dashboard (Shiny/Dash) for exploratory analysis
- Docker/Singularity containerisation for HPC/cloud deployment

---

## References

1. Feldgarden M et al. (2021) AMRFinderPlus. *Scientific Reports* 11:12115
2. Katoh K, Standley DM (2013) MAFFT. *Mol Biol Evol* 30:772-780
3. Seemann T (2014) Prokka. *Bioinformatics* 30:2068-2069
4. Unemo M, Shafer WM (2014) Antimicrobial resistance in *Neisseria gonorrhoeae*. *Clin Microbiol Rev* 27:587-613
5. NCBI Reference: GCF_000006845.1 - *Neisseria gonorrhoeae* FA 1909, complete genome
