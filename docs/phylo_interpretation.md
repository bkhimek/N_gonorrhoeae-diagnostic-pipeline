# Phylogenetic Analysis — *N. gonorrhoeae* AMR Interpretation

**Pipeline:** Parsnp v2.1.5 · FA1909 reference · 284 genomes (283 strains + reference)  
**AMR annotation:** AMRFinderPlus v4.0.23 · 277 strains with resistance data  
**Date:** June 2026

---

## 1. Overview

A core-genome SNP phylogeny was reconstructed from 283 *N. gonorrhoeae* assemblies
using Parsnp v2.1.5 with FA1909 as the reference strain. Tip points are coloured
by AMR gene burden (number of unique resistance elements per strain, range 0–20).
The interactive tree is available as `phylo_amr_tree.html`.

---

## 2. Phylogenetic Structure

The tree reveals clear clade structure inconsistent with a star phylogeny, indicating
that the 283 strains represent genuine evolutionary divergence rather than a single
outbreak clone. Several well-supported lineages are visible:

- **Upper clade (large, diverse):** The most speciose lineage, containing the majority
  of strains. Internal branching is shallow, consistent with recent diversification.
- **Mid-tree lineages:** Two to three intermediate-sized clades with moderate branch
  lengths, suggesting distinct but related genetic backgrounds.
- **Basal clade (bottom):** A well-separated lineage with notably lower AMR burden
  (predominantly yellow/cream tips), possibly representing an older or less clinically
  prevalent genomic background.
- **Long-branch outliers:** A small number of strains with branch lengths substantially
  longer than their neighbours may represent more divergent lineages, rarer sequence
  types, or assembly artefacts worth flagging for manual review.

---

## 3. AMR Burden Distribution

### 3.1 Clustering of high-burden strains

High AMR burden strains (≥15 unique resistance elements, dark red) are **not randomly
distributed** across the tree. They cluster preferentially within specific subclades,
most prominently in the upper portion of the tree. This pattern is consistent with
**clonal expansion of multi-drug resistant (MDR) lineages** — a small number of
successful genotypes have spread while accumulating resistance determinants, rather
than resistance arising independently many times.

### 3.2 Gradient lineages

Several clades show a visible gradient from low (yellow) to high (red) burden along
terminal branches. This likely reflects **stepwise acquisition of resistance** within
an evolving lineage — strains within the same genomic background have progressively
accumulated resistance elements over time, potentially driven by antibiotic selection
pressure.

### 3.3 Low-burden basal clade

The basal clade contains predominantly low-burden strains (0–5 elements). If these
strains represent an earlier-diverging lineage, this may indicate that the MDR
phenotype is a more recent evolutionary acquisition, emerging after the major
lineage split.

### 3.4 Highly resistant outliers

A cluster of strains with burden ≥18–20 elements (deepest red) appears as a tight
sub-clade in the mid-upper tree. These likely represent **extensively drug-resistant
(XDR)** candidates carrying simultaneous beta-lactam, fluoroquinolone, tetracycline,
sulfonamide, and rifamycin resistance determinants.

---

## 4. Drug-Class Highlights

From the AMR dashboard (drug class summary tab):

| Drug class | Strains affected | % |
|---|---|---|
| Beta-lactam | 264 | 95% |
| Tetracycline | 211 | 76% |
| Sulfonamide | 201 | 73% |
| Quinolone | 171 | 62% |
| Rifamycin | 86 | 31% |
| Beta-lactam/Macrolide/Tetracycline | 116 | 42% |

Beta-lactam resistance is near-universal, consistent with the very high prevalence
of *penA* alleles (penA_F504L and penA_A510V each at ~80%). Quinolone resistance
(62%) is driven primarily by *gyrA_S91F* and *folP_R228S* point mutations. The
multi-class BETA-LACTAM/MACROLIDE/TETRACYCLINE category (42% of strains) likely
reflects *mtrR_A39T*, an efflux pump promoter mutation conferring broad resistance.

---

## 5. Limitations

- **Core genome only:** Parsnp aligns only the conserved core genome. Accessory
  resistance elements (plasmid-borne genes, genomic islands) present in some strains
  may be underrepresented in the SNP alignment.
- **Reference bias:** Using FA1909 as reference may introduce bias if some strains
  are highly divergent from this background.
- **No bootstrap support:** The Parsnp/FastTree tree is unrooted/midpoint-rooted
  without formal bootstrap resampling. Clade support should be interpreted
  cautiously without further validation (e.g., IQ-TREE with bootstrap).
- **AMR burden vs. phenotype:** Gene burden is a proxy for resistance; phenotypic
  MIC data would be required to confirm clinical resistance levels.

---

## 6. Next Steps

1. **Bootstrap the tree** using IQ-TREE2 for clade support values.
2. **Metadata overlay** — if collection date/country/source metadata are available,
   overlay on tree tips to assess spatiotemporal spread of MDR clones.
3. **Promoter-region variants** — extend AMR analysis to *mtrR*, *penB*, and *ponA*
   promoter regions to capture regulatory resistance mechanisms not detected by
   AMRFinderPlus.
4. **Plasmid analysis** — screen assemblies for known *N. gonorrhoeae* plasmids
   (e.g., *β*-lactamase plasmids pβ1, pβ2) to complement chromosomal AMR findings.

---

*Generated with the N. gonorrhoeae AMR diagnostic pipeline.*  
*Repository: https://github.com/bkhimek/N_gonorrhoeae-diagnostic-pipeline*
