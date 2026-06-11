# Lessons Learned: N. gonorrhoeae Diagnostic Pipeline

## Overview

This document captures the practical challenges, analytical pitfalls, and hard-won insights from building and running the *N. gonorrhoeae* diagnostic target discovery pipeline. It is intended as a reference for running similar pipelines on other bacterial species.

---

## 1. The Biggest Analytical Mistake: Per-Hit vs. Union Coverage

**What happened:** The initial exclusivity analysis reported that 16 genes passed, with per-hit query coverage (qcovs) of 51–77%. This looked promising — it seemed like each gene had large regions not covered by any BLAST hit.

**Why it was wrong:** qcovs is per-hit. Different organisms hit different parts of the gene. A hit from *N. sicca* might cover positions 200–600; a hit from *N. lactamica* covers positions 400–900. Neither hit alone covers the full gene, but together they cover nearly all of it. Reading per-hit qcovs as "the uncovered flanks" is a fundamental mistake.

**The fix:** Merge all BLAST hit intervals per gene (union-of-hits), then find gaps in the merged set. This is the only correct way to identify genuinely uncovered positions.

**Lesson:** Always perform interval merging when assessing coverage from multiple overlapping BLAST hits. Per-hit statistics are misleading in exactly this scenario.

```python
# Correct approach: merge intervals before finding gaps
intervals = sorted(covered[gene])
merged = []
for start, end in intervals:
    if merged and start <= merged[-1][1]:
        merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    else:
        merged.append([start, end])
```

---

## 2. The N. gonorrhoeae / N. meningitidis Problem

**The challenge:** *N. gonorrhoeae* and *N. meningitidis* diverged relatively recently and share ~80–85% genome-wide nucleotide identity. This is the fundamental biological obstacle to gonorrhoea-specific diagnostics.

**What it means in practice:** Most gonococcal genes have meningococcal orthologues with BLAST hits at e-value < 1e-5. After interval merging, 12 of 16 candidate PASS genes had zero uncovered sequence. The species simply does not have many genes that are absent or highly diverged in *N. meningitidis*.

**Lesson for other species:** Before running this pipeline, assess phylogenetic distance to the nearest non-target relative. If that relative is within the same genus, expect few or no unique regions among housekeeping genes. You may need to look at accessory genome elements (pathogenicity islands, unique virulence factors) rather than core genes.

---

## 3. Excluding Target Taxids from BLAST: Use `-negative_taxidlist`

**The problem:** BLAST's `-negative_taxids` parameter accepts a comma-separated list, which breaks when you have more than ~100 taxids. *N. gonorrhoeae* has 485 (species node) plus 89+ strain-level child taxids in NCBI Taxonomy.

**The solution:** Use `-negative_taxidlist` with a plain text file (one taxid per line). This works reliably for any number of taxids.

**How to get all child taxids:**
```bash
# Using NCBI E-utilities
efetch -db taxonomy -id 485 -format xml | \
    xtract -pattern Taxon -element TaxId > gonorrhoeae_all_taxids.txt
```

**Lesson:** Always use `-negative_taxidlist` (file) rather than `-negative_taxids` (inline) for any species with more than a handful of strain entries in NCBI Taxonomy. Most pathogens with significant research interest will have dozens to hundreds of strain-level taxids.

---

## 4. Pan-genome Size and Core Gene Threshold

**Observation:** With 283 genomes, the *N. gonorrhoeae* pan-genome is large. The core genome (≥95% prevalence) contains far fewer genes than the full pan-genome. This is by design — you want genes present in essentially all strains.

**Tradeoff:** A higher threshold (e.g., 99%) gives a smaller, safer target list; a lower threshold (e.g., 90%) gives more candidates but risks missing the target in some strains. For a diagnostic assay, ≥95% is a reasonable starting point, and any selected target should be verified against strain collections before clinical use.

**Lesson:** Report the exact threshold used and the number of core genes it yields. Different Panaroo runs with different thresholds can give very different candidate lists.

---

## 5. BLAST Database Choice Matters

**Used:** `nt_prok` (prokaryotic sequences only from NCBI nt)

**Why not full `nt`:** The full nt database is enormous (~500 GB+). For exclusivity screening of bacterial targets, nt_prok is sufficient and dramatically faster.

**Caveat:** nt_prok excludes eukaryotic sequences. If the intended clinical application involves specimens that might contain human DNA or fungal sequences, a full nt BLAST (or a separate specificity check against human genome) may be warranted.

**Lesson:** Match the BLAST database to your biological question. For species-level bacterial exclusivity: nt_prok. For clinical specificity including host background: full nt or add a separate human BLAST check.

---

## 6. Panaroo vs. Roary

**Why Panaroo was chosen:** Roary is faster but known to overestimate the core genome (it merges gene families too aggressively). Panaroo uses a graph-based approach that is more conservative and produces fewer spurious core gene calls. For diagnostic work where false positives (thinking a gene is in all strains when it isn't) are costly, Panaroo is the safer choice.

**Lesson:** For clinical/diagnostic pipeline development, err on the side of conservative pan-genome tools. A missed candidate is recoverable; a false positive wastes wet-lab resources.

---

## 7. File Paths and Docker Volumes

**Pain point:** The Nextflow pipeline runs tools inside Docker containers. File paths inside the container differ from host paths. This caused silent failures when BLAST databases or genome files were not accessible inside the container.

**Fix:** Explicitly mount all required directories in `nextflow.config`:
```groovy
docker {
    enabled = true
    runOptions = '-v /home/himek/blast_db:/blast_db -v /home/himek/genomes:/genomes'
}
```

**Lesson:** When setting up a Dockerised Nextflow pipeline, list every external directory the pipeline needs and add it as a volume mount. Test with a single genome before scaling.

---

## 8. Representative Sequence Selection for Primer Design

**Issue:** The per-gene FASTA files from `pass_fastas/pass/` contain one sequence per input genome (283 sequences per gene). Primer3 needs a single template.

**Approach taken:** Used the first sequence in the FASTA as a representative. This is adequate for an initial screen but not ideal.

**Better approach for production:** Compute a multiple sequence alignment of all 283 sequences for the unique region, identify conserved positions, and design primers only to fully conserved positions. SNPs in the primer binding site can cause allele dropout in PCR.

**Lesson:** For any diagnostic assay, always check primer binding sites against the full strain diversity of your target organism before going to the lab.

---

## 9. Interpreting "PASS" from the Pipeline

**What PASS means:** The gene had at least one BLAST hit with qcovs < some threshold (or no hit at all to non-gonococci), suggesting some part of it may be unique.

**What PASS does NOT mean:** The gene has a unique region large enough for PCR. As shown here, 14 of 16 PASS genes had zero usable unique sequence after interval merging.

**Lesson:** The pipeline's PASS/FAIL is a first-pass filter. The real analysis is the interval-merging step. Don't skip it.

---

## 10. Challenging Target = Good Portfolio Piece

**Reframe:** The fact that only 2 genes (polA, norR) emerged as viable targets out of 283 core genes is not a failure — it is the scientifically correct answer to a genuinely hard problem. It reflects the real biological constraint that gonococci and meningococci are nearly identical genomically.

**Portfolio value:** A case study that shows you identified a subtle analytical error (per-hit vs. union coverage), corrected it, and arrived at a well-justified, conservative conclusion is more impressive than one where everything works cleanly from the start.

---

## Summary Checklist for Running This Pipeline on a New Species

- [ ] Confirm phylogenetic distance to nearest non-target relative
- [ ] Retrieve all strain-level taxids for the target species
- [ ] Use `-negative_taxidlist` (not `-negative_taxids`) in BLAST config
- [ ] Run Panaroo at ≥95% core threshold
- [ ] After pipeline: run interval-merging analysis (not per-hit qcovs)
- [ ] Filter for unique regions ≥150 bp
- [ ] Use multi-sequence alignment for primer binding site conservancy check
- [ ] Validate top candidates with NCBI Primer-BLAST before wet lab
