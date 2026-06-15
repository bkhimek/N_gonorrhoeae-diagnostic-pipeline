# Promoter-Region Variant Analysis — Interpretation

*N. gonorrhoeae* diagnostic pipeline · 283 strains · FA1090 reference (NC_002946.2)

---

## 1. Overview

Resistance to multiple antibiotic classes in *Neisseria gonorrhoeae* is mediated not only by
acquired resistance genes but also by regulatory mutations that upregulate intrinsic efflux
pumps. The most clinically important regulatory locus is the **mtrR/mtrCDE intergenic region**,
where a 13-bp inverted repeat (IR; 5′-TAATTTTTATTTT-3′) acts as an operator. Deletion of this
IR derepresses the MtrCDE efflux pump, increasing resistance to penicillin, tetracycline,
and macrolides (Hagman & Shafer, 1995; PMID 7506352).

This analysis extracted the 400 bp upstream of *mtrR* and 300 bp upstream of *porB* from all
283 assemblies, aligned each to the FA1090 reference region using MAFFT, and quantified
single-nucleotide polymorphisms (SNPs) and insertion/deletion events (indels).

---

## 2. Reference Genome Note

FA1090 (GCF_000006845.1), the standard *N. gonorrhoeae* reference, **does not contain** the
classic 13-bp IR sequence in its *mtrR* upstream region. Consequently, clinical strains that
*retain* the intact IR appear as insertions relative to FA1090, rather than the expected
deletions. All variant calls in this analysis are therefore expressed as differences from the
FA1090 promoter state.

---

## 3. mtrR Upstream Region

### 3.1 BLAST coverage

All 283 assemblies produced a BLAST hit against the FA1090 *mtrR* upstream reference (100%
detection rate, average identity ≥ 99%).

### 3.2 Indel distribution

| Indels vs FA1090 | Strains | % of cohort | Interpretation |
|:---:|---:|---:|---|
| 0 | 1 | 0.4 | FA1090 reference itself (identical) |
| 1 | 157 | 55.5 | Minor structural variant (most common) |
| 2 | 112 | 39.6 | Two-position structural change |
| 3 | 4 | 1.4 | Small insertion/deletion |
| 5 | 1 | 0.4 | Small insertion/deletion |
| 8 | 2 | 0.7 | Moderate structural variant |
| 9 | 1 | 0.4 | Moderate structural variant |
| 27 | 1 | 0.4 | ⚠ Possible assembly issue |
| 189 | 4 | 1.4 | ⚠ Likely assembly artefact |

The dominant pattern (1–2 indels, 269/283 strains = 95.1%) is consistent with normal
sequence polymorphism in an AT-rich regulatory region. Whether these indels correspond to
the presence or absence of the 13-bp IR cannot be resolved with FA1090 as the sole
reference, because FA1090 itself already lacks the IR.

### 3.3 SNP distribution

The majority of strains carry 0 or 1 SNP relative to FA1090 in the *mtrR* upstream region
(128 and 141 strains respectively), indicating that point mutations in this promoter are
common but modest. A small number of outlier strains (≤ 7) carry 47–77 SNPs, consistent
with the same assembly-quality issues that inflate their indel counts.

### 3.4 Flagged strains

Five strains exceeded the quality threshold of 20 indels:

| Strain | mtrR indels | AMR elements | Assessment |
|--------|:-----------:|:------------:|------------|
| GCF_040373105.1 | 189 | 8 | Localised assembly artefact |
| GCF_900087735.2 | 189 | 8 | Localised assembly artefact |
| GCF_900637245.1 | 189 | 8 | Localised assembly artefact |
| GCF_978047055.1 | 189 | 8 | Localised assembly artefact |
| GCF_978047425.1 | 27 | 8 | Possible localised assembly artefact |

All five carry 8 AMR elements — below the cohort median (12) but within the expected range
(0–23). Their overall genome quality appears intact; the inflated indel count is confined to
the *mtrR* locus and likely reflects difficulty in assembling this AT-rich, repetitive region.
These strains are **retained in the main AMR analysis** but should be interpreted with caution
for promoter-level calls.

---

## 4. porB Upstream Region

| Indels vs FA1090 | Strains | % of cohort |
|:---:|---:|---:|
| 0 | 276 | 97.5 |
| 1 | 4 | 1.4 |
| 160 | 3 | 1.1 |

The *porB* upstream region is highly conserved across this cohort: 97.5% of strains show
no indels relative to FA1090. The three strains with 160 indels are the same assembly-quality
outliers identified in the *mtrR* analysis. No specific known *porB* promoter variants
were called in this dataset.

---

## 5. Limitations

1. **Reference choice.** Using FA1090 (which lacks the IR) as the reference means IR-intact
   (susceptible) and IR-deleted (resistant) strains cannot be distinguished by indel
   direction alone. A susceptible-strain reference (e.g. FA19, WHO reference F) would allow
   direct IR deletion calling.

2. **Alignment sensitivity.** MAFFT pairwise alignment of short AT-rich regions can
   misplace gap positions, affecting indel counts in edge cases.

3. **Indel counting method.** The reported "indels" are alignment gap-columns, not the
   number of distinct insertion/deletion events; a single 13-bp deletion counts as 13 indels.

4. **No phenotypic data.** Without MIC measurements, resistance/susceptibility phenotypes
   cannot be confirmed from promoter genotype alone.

---

## 6. Next Steps

- **Use a susceptible reference.** Download a susceptible FA19-type strain and use its
  *mtrR* upstream region as the reference to enable direct IR deletion detection.
- **Validate outlier strains.** If short-read Illumina data are available for the five
  flagged strains, map reads to confirm the assembly call in this region.
- **Correlate with MIC data.** Cross-reference promoter indel status with phenotypic MIC
  data to validate the resistance prediction.
- **Integrate with AMR gene data.** Strains with *mtrR* promoter variants and concurrent
  *penA*, *gyrA*, and *parC* mutations represent high-priority MDR candidates.

---

## 7. References

- Hagman KE, Shafer WM (1995). Transcriptional mechanism of the *mtrCDE* efflux pump.
  *J Bacteriol* 177(14):4162–4165. PMID 7506352.
- Unemo M, Shafer WM (2014). Antimicrobial resistance in *Neisseria gonorrhoeae*.
  *Clin Microbiol Rev* 27(3):587–613. PMID 24982323.
