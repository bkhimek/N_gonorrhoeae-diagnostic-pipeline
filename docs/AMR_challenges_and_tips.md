# AMR Module — Challenges, Fixes, and Reproducibility Tips

*N. gonorrhoeae* Diagnostic Pipeline | June 2026

---

## Introduction

This document records every significant issue encountered during the build of the AMR resistance gene detection module, the fix applied in each case, and generalised tips for reproducing or adapting this pipeline for other bacterial species. Issues are grouped by category.

---

## 1. AMRFinderPlus Configuration

### 1.1 Wrong organism name
**Issue:** `--organism Neisseria` raised an error — AMRFinderPlus requires the full binomial.  
**Fix:** Use `--organism Neisseria_gonorrhoeae` (or `Neisseria_meningitidis`). Run `amrfinder --list_organisms` to see valid names.  
**Tip for other species:** Always check the supported organism list first. If your species is absent, run without `--organism` — you still get gene-level detection, just no organism-specific point mutation rules.

### 1.2 Database cannot be updated to a custom path
**Issue:** `amrfinder -u --database /custom/path` fails — the update option only works on the default database directory.  
**Fix:** Run `amrfinder -u` without `--database`. The database installs to `~/miniconda3/share/amrfinderplus/data/latest`. Reference this path in `nextflow.config` using `System.getenv('HOME')` instead of `$HOME` (Nextflow does not expand shell variables in config files directly).  
**Tip:** Pin the database version by referencing the dated subfolder (e.g., `data/2025-07-16.1`) in production runs so results are reproducible even after database updates.

### 1.3 GFF format incompatibility with Prokka
**Issue:** Running AMRFinderPlus in full mode (`--protein --gff --nucleotide`) fails with `GFF file mismatch: Protein FASTA id not in GFF file`. Prokka GFF3 uses `ID=`/`locus_tag=` attributes; AMRFinderPlus expects `protein_id=` (GenBank-style).  
**Fix:** Run AMRFinderPlus in **protein-only mode** (`--protein` only, no `--gff` or `--nucleotide`). This still detects resistance genes and protein-level point mutations (POINTP method). Clinical point mutations are handled by the custom alignment pipeline.  
**Note:** Stripping the `##FASTA` block from Prokka GFF (`sed '/^##FASTA/,$d'`) resolves the line-count error but not the attribute mismatch — protein-only mode is the correct long-term solution unless you reformat the GFF.  
**Tip for other species:** If using NCBI-annotated GFF (not Prokka), the full AMRFinderPlus mode should work without modification.

---

## 2. BLAST Gene Extraction

### 2.1 BLAST database names with dots cause silent failures
**Issue:** Creating a BLAST database with a dot in the filename (e.g., `/tmp/blast_GCF_000006845.1`) combined with `-quiet` redirected to `/dev/null` caused silent failure — blastp returned no hits even though the query gene was present.  
**Fix:** Replace dots with underscores in database names (`${sid//./_}` in bash). Remove `-quiet` during debugging. The Nextflow module uses `strain_db` as a safe generic name within the work directory.  
**Tip:** Never use `-quiet 2>/dev/null` during development — errors are hidden and debugging becomes very difficult.

### 2.2 `samtools faidx` fails silently on unindexed files
**Issue:** Extracting a BLAST hit sequence with `samtools faidx $faa "$hit"` produces no output if the `.fai` index does not exist.  
**Fix:** Run `samtools faidx $faa` to create the index before attempting random-access extraction.  
**Tip:** In Nextflow, include `samtools faidx` as part of the process script rather than assuming pre-indexed input.

### 2.3 BLAST `-max_target_seqs 1` warning
**Issue:** `Warning: [blastp] Examining 5 or more matches is recommended` when using `-max_target_seqs 1`.  
**Fix:** Use `-max_target_seqs 5` and filter the output with `awk` to take the top hit passing identity/coverage thresholds. This also avoids the known BLAST behaviour where `-max_target_seqs 1` does not guarantee the best hit is returned.

---

## 3. MAFFT

### 3.1 MAFFT not installed in base conda environment
**Issue:** `mafft: command not found` — the base conda environment did not include MAFFT.  
**Fix:** `conda install -y -c bioconda mafft`  
**Tip for reproducibility:** Create a dedicated conda environment (e.g., `n_gonorrhoeae_amr`) and capture it with `conda env export > environment.yml`. This ensures all tools are versioned and the environment can be recreated on any system.

---

## 4. Nextflow DSL2 Configuration

### 4.1 `$HOME` not resolved in `nextflow.config`
**Issue:** Using `"$HOME/path/to/db"` in `nextflow.config` raises `HOME is not defined (hint: use env('...') to access environment variable)` in Nextflow 26.x.  
**Fix:** Use `System.getenv('HOME') + '/path/to/db'` in the config file.

### 4.2 `fromPath` with glob requires `type: 'dir'`
**Issue:** `Channel.fromPath("prokka_output/*/")` returns an empty channel in Nextflow 26 when the glob targets directories.  
**Fix:** Add `type: 'dir'`: `Channel.fromPath("prokka_output/*/", type: 'dir')`.  
**Tip:** Always verify channel content during development with `.view()` before connecting to processes.

### 4.3 `projectDir` resolves to the script's location, not the project root
**Issue:** When running `nextflow run /tmp/test.nf`, `projectDir` is `/tmp/`, so all `${projectDir}/reference/...` paths point to non-existent locations. Processes are silently dropped rather than erroring.  
**Fix:** Always run `nextflow run` from the project root directory. For testing, use a test script that imports from absolute paths, or pass `--reference_dir` as a command-line parameter.  
**Tip:** This is one of the most common silent failure modes in Nextflow — if processes show `-` (not started) instead of errors, check that all input paths actually exist.

### 4.4 Collected file lists cannot be passed as directory arguments
**Issue:** When `CALL_MUTATIONS.out` files are collected and passed to `RESISTANCE_PROFILE`, Nextflow expands them as a space-separated list of filenames. A Python script expecting a directory path receives multiple positional arguments instead.  
**Fix:** In the process script, move collected files into a subdirectory before calling the Python script:
```bash
mkdir -p mut_dir
for f in *_mutations.tsv; do mv "$f" mut_dir/; done
python3 summarize.py --mutations-dir mut_dir
```

### 4.5 Complex channel operations in subworkflow
**Issue:** Attempting to write files inside a Nextflow `.map{}` closure (side effects) and using multi-step `.flatMap`/`.collect` chains caused processes to silently not run.  
**Fix:** Move file-writing logic into a dedicated `COLLECT_GENES` process. Keep channel operations purely transformational (no side effects in closures). When in doubt, break complex channel transformations into intermediate processes.

---

## 5. Python Scripts

### 5.1 Long heredocs fail in WSL terminal
**Issue:** Pasting a heredoc (`cat > file << 'EOF' ... EOF`) longer than ~50 lines into a WSL terminal results in partial or no file creation — the terminal buffer is exceeded.  
**Fix:** Split the file into chunks of ~20–25 lines, using `cat >` for the first chunk and `cat >>` for subsequent chunks.  
**Tip:** For files >50 lines, write them on the Windows side (e.g., VS Code) and access via `/mnt/c/...` from WSL, or use a Python script that writes the content programmatically.

### 5.2 `awk` cannot sum Python `True`/`False` strings
**Issue:** `awk '{sum += $16}'` on a column containing `"True"`/`"False"` strings produces 0 — awk converts non-numeric strings to 0.  
**Fix:** Use string comparison: `if($16=="True") count++`. Or use Python's `csv.DictReader` to read by column name, which is more robust.

### 5.3 Reference sequence not found in alignment
**Issue:** `parse_mutations.py` failed with `ERROR: reference not found` because the reference FASTA header was renamed to `>gyrA_FA1909` (not `>GCF_000006845.1`).  
**Fix:** Pass `--reference FA1909` (the string present in the header) rather than the full GCF accession. Consistent header naming between the extraction and calling steps is essential.

---

## 6. General Reproducibility Tips for Other Bacteria

### Choosing a reference strain
- Select a fully closed, well-annotated reference genome (ideally a complete chromosome, not a draft assembly).
- Confirm the genome is widely used in the literature so your mutation numbering matches published studies.
- For *N. gonorrhoeae*, FA-1909 (GCF_000006845.1) is the standard. For *N. meningitidis*, use H44/76 (GCF_000008805.1).

### Identifying target gene locus tags
- Prokka assigns arbitrary locus tags — do not assume gene names are identical across strains or annotation runs.
- Use the `.gbk` file's `/gene=` attribute (not just the product name) to reliably identify target genes: `grep -A2 'gene="gyrA"' *.gbk`
- Cross-check locus tag → protein length against known values for your organism.

### Setting BLAST thresholds
- The 80% identity / 80% query coverage defaults used here are appropriate for well-conserved housekeeping genes within a single species.
- For divergent genes or multi-species datasets, lower thresholds (e.g., 60% identity) may be needed, but increase the risk of paralogue hits.
- Always inspect the raw BLAST output for at least a subset of strains to confirm the correct gene is being retrieved.

### QRDR window definitions
- QRDR positions are species- and gene-specific. The windows used here (GyrA aa 67–106, ParC aa 57–96) are validated for *Neisseria*.
- For other species (e.g., *Campylobacter*, *Salmonella*), consult the literature for the correct QRDR boundaries before adapting `parse_mutations.py`.

### Validating the pipeline
- Include at least two control strains with known resistance profiles (e.g., a susceptible WHO F strain and a multi-drug resistant WHO K strain).
- Confirm that `blaTEM-1` strains are detected by AMRFinderPlus and that known ciprofloxacin-resistant strains show GyrA S91F.

### Version pinning
- Pin tool versions in `nextflow.config` or a `conda environment.yml`.
- Record the AMRFinder database version (dated subfolder) in your methods — the database is updated frequently and results can change.
- Use `nextflow -resume` for re-runs to avoid recomputing completed steps when only one module changes.
