# Case Study: Developing Species-Specific PCR Diagnostic Targets for *Neisseria gonorrhoeae*

**Domain:** Bioinformatics | Molecular Diagnostics | Bacterial Genomics  
**Tools:** Nextflow DSL2, Panaroo, BLAST+, Primer3, Python  
**Data:** 283 complete *N. gonorrhoeae* RefSeq genomes; NCBI nt_prok database

---

## Background

*Neisseria gonorrhoeae* (the gonococcus) is the causative agent of gonorrhoea, a globally prevalent sexually transmitted infection with rising antimicrobial resistance. Rapid, accurate molecular diagnosis is critical for patient management and surveillance. The challenge in designing species-specific PCR assays is extreme: *N. gonorrhoeae* is phylogenetically almost indistinguishable from its commensal relative *N. meningitidis*, which colonises the nasopharynx of healthy individuals. Any diagnostic target must be present in essentially all gonococci (inclusivity) while producing no signal from meningococci or any other organism (exclusivity).

This case study describes a fully automated bioinformatics pipeline to identify and validate candidate PCR target genes, culminating in the extraction of primer-ready sequences for the top candidates.

---

## Objectives

1. Identify genes present in ≥95% of *N. gonorrhoeae* genomes (core genome, inclusivity)
2. Identify which of those genes have regions with no significant BLAST hits to non-gonococcal organisms (exclusivity)
3. Determine which exclusive regions are large enough (≥150 bp) to accommodate a PCR primer pair
4. Extract the sequences of those regions and design candidate primers with Primer3

---

## Pipeline Overview

```
283 N. gonorrhoeae genomes (RefSeq)
        │
        ▼
   [Prokka] — Gene annotation
        │
        ▼
   [Panaroo] — Pan-genome analysis
        │  Core genome genes (≥95% prevalence)
        ▼
   [BLAST exclusivity] — vs. nt_prok, excluding all 90 N. gonorrhoeae taxids
        │  (-negative_taxidlist, e-value 1e-5)
        ▼
   [Coverage analysis] — Union-of-hits interval merging per gene
        │
        ▼
   Candidate unique regions → Primer3
```

### Key pipeline parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Core genome threshold | ≥95% genomes | Ensures broad strain coverage |
| BLAST e-value cutoff | 1e-5 | Standard for homology detection |
| Excluded taxids | 90 (taxid 485 + all strain children) | Remove all gonococci from BLAST db |
| Min unique region size | ≥150 bp | Minimum for functional primer pair |
| Primer Tm range | 57–63°C (opt 60°C) | Standard PCR conditions |
| Amplicon size range | 150–400 bp | Suitable for conventional PCR |

---

## The Core Biological Challenge

The most important finding emerged early: *N. gonorrhoeae* and *N. meningitidis* share ~80–85% average nucleotide identity across their genomes. This means that most gonococcal genes have meningococcal orthologues, and those orthologues have high enough similarity to produce BLAST hits under standard thresholds.

A naive reading of the BLAST results was misleading. Sixteen genes passed an initial exclusivity filter (no single BLAST hit covered >77% of the gene). However, analysis of per-hit query coverage (qcovs 51–77%) initially suggested large uncovered flanks — which would be ideal primer targets.

**The critical insight:** different organisms hit *different parts* of the same gene. When all BLAST hits are merged into a union of covered intervals, most genes are fully covered. The per-hit flanks were an artefact of the analysis method, not genuine unique sequence.

---

## Interval Merging Analysis

To find genuinely uncovered regions, BLAST coordinates (qstart, qend) were extracted from `blast_results.tsv` for all 16 candidate PASS genes. For each gene, all hit intervals were merged using the following algorithm:

```python
intervals = sorted(covered[gene])
merged = []
for start, end in intervals:
    if merged and start <= merged[-1][1]:
        merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    else:
        merged.append([start, end])
# gaps between merged intervals = truly uncovered regions
```

### Results

| Gene | Gene Length | Largest Uncovered Region | Size | PCR Suitability |
|------|-------------|--------------------------|------|-----------------|
| **polA** | 2823 bp | pos 1–910 (5′ end) | **910 bp** | ✓ Primary target |
| **norR** | 1727 bp | pos 1015–1727 (3′ end) | **713 bp** | ✓ Primary target |
| **pyrI** | 525 bp | pos 1–161 (5′ end) | 161 bp | Marginal |
| bcp | 435 bp | pos 1–1 | 1 bp | Not suitable |
| 12 others | — | 0 bp | 0 bp | Not suitable |

Only 2 genes (polA, norR) have unique regions of sufficient size for reliable PCR primer design. This result is a direct consequence of the extreme phylogenetic proximity of *N. gonorrhoeae* to *N. meningitidis*.

---

## Candidate Target Genes

### polA (DNA polymerase I) — Primary Target

- **Unique region:** 5′ first 910 bp of a 2823 bp gene
- **Biological context:** DNA polymerase I, involved in replication and repair. The N-terminal domain diverges between gonococci and meningococci despite the central polymerase domain being conserved
- **PCR design:** Ample sequence for multiple primer pairs; 150–400 bp amplicons easily achievable
- **Confidence:** High — large unique window, core gene present in all 283 genomes

### norR (nitric oxide-responsive transcriptional regulator) — Primary Target

- **Unique region:** 3′ terminal 713 bp (pos 1015–1727)
- **Biological context:** Transcriptional regulator of the nitric oxide stress response (norVW operon). Divergence in the C-terminal effector domain likely reflects adaptation of gonococci to the specific oxidative environment of the urogenital tract
- **PCR design:** Good size; 3′ location means primers can span the unique region cleanly
- **Confidence:** High

### pyrI (carbamoyl-phosphate synthase small subunit) — Marginal

- **Unique region:** 5′ first 161 bp of a 525 bp gene
- **Biological context:** Pyrimidine biosynthesis. Short unique window limits primer placement flexibility
- **PCR design:** Borderline — 161 bp leaves little room; a single primer pair at the extreme 5′ end is theoretically possible but robustness is uncertain
- **Confidence:** Low — not recommended for primary assay development

---

## Primer Design (Primer3)

Primer3 was run on each unique region with the following constraints:

```
PRIMER_OPT_SIZE=20, PRIMER_MIN_SIZE=18, PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0, PRIMER_MIN_TM=57.0, PRIMER_MAX_TM=63.0
PRIMER_MIN_GC=40.0, PRIMER_MAX_GC=65.0
PRIMER_PRODUCT_SIZE_RANGE=150-400
PRIMER_NUM_RETURN=3
```

The sequences used as input are provided in `unique_regions_for_primers.fasta`.

**Recommended next step:** Submit polA and norR primer candidates to NCBI Primer-BLAST with organism filter set to *N. gonorrhoeae* (taxid 485) to confirm in-silico specificity before wet-lab validation.

---

## Limitations and Future Work

**N. meningitidis specificity gap:** The exclusivity BLAST search excluded all *N. gonorrhoeae* taxids but included *N. meningitidis*. A dedicated pairwise alignment between gonococcal and meningococcal polA/norR sequences would quantify the exact divergence in the unique regions and strengthen confidence in assay specificity.

**Single representative sequence:** Primer design was performed on one representative genome per gene. A multi-sequence alignment of all 283 gonococci for the unique region would identify SNPs that could destabilise primer binding across strain diversity.

**Wet-lab validation required:** Computational predictions must be validated by:
1. In-silico PCR (e-PCR / Primer-BLAST)
2. PCR on a panel of *N. gonorrhoeae* strains (sensitivity)
3. PCR on *N. meningitidis*, other *Neisseria* spp., and common urogenital flora (specificity)

**AMR marker integration:** Resistance-associated mutations in penA, mtrR, porB, gyrA, and parC could be incorporated as secondary targets in a multiplex format to provide simultaneous diagnosis and resistance profiling.

---

## Conclusion

Of 283 core *N. gonorrhoeae* genes assessed, only two — **polA** and **norR** — contain regions genuinely absent from all non-gonococcal organisms in the nt_prok database. The extreme genomic similarity between *N. gonorrhoeae* and *N. meningitidis* means that most candidate genes are disqualified once all BLAST hits across organisms are merged. This finding underscores the difficulty of designing gonorrhoea-specific diagnostics and highlights polA and norR as the most promising computational targets for wet-lab follow-up.

---

## Files

| File | Description |
|------|-------------|
| `unique_regions_for_primers.fasta` | Extracted unique region sequences for polA, norR, pyrI |
| `results/exclusivity/blast_results.tsv` | Full BLAST exclusivity results (on himek@MSI) |
| `results/exclusivity/exclusivity.tsv` | Pipeline PASS/FAIL calls per gene |
| `results/pass_fastas/pass/` | Per-gene FASTA files for all PASS genes |
| `unique_regions.py` | Interval-merging coverage analysis script |
