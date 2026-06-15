#!/usr/bin/env python3
"""
detect_promoter_variants.py — Detect resistance-associated promoter variants.

Uses BLAST to locate promoter regions in each assembly, BioPython for
sequence extraction, MAFFT for alignment, then calls known variants.

Key variant:
  mtrR promoter — 13-bp inverted repeat (IR) deletion → upregulates mtrCDE.
  Susceptible strains carry the IR (TAATTTTTATTTT); resistant strains lack it.

  IMPORTANT: FA1090 (reference) appears to lack the classic 13-bp IR.
  We therefore search for the IR directly in each strain's promoter sequence
  rather than using absence-vs-reference logic.

Usage (single strain):
    python detect_promoter_variants.py \
        --assembly  genomes/GCF_001234.fna \
        --promoters ~/n_gonorrhoeae_project/reference/promoters/ \
        --out       results/promoter/GCF_001234_promoter.tsv
"""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# ── Known signatures ──────────────────────────────────────────────────────────
# 13-bp IR present in susceptible strains (Hagman & Shafer 1995, PMID 7506352).
# Deletion of this IR upregulates mtrCDE efflux pump → resistance.
# Search both forward and reverse complement.
MTRR_IR    = "TAATTTTTATTTT"
MTRR_IR_RC = str(Seq(MTRR_IR).reverse_complement())  # AAAATAAAATTA


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout.strip(), r.stderr.strip()


def load_fasta_index(fna):
    return SeqIO.to_dict(SeqIO.parse(fna, "fasta"))


def match_contig(records, blast_id):
    if blast_id in records:
        return records[blast_id]
    for rid, rec in records.items():
        if blast_id in rid or rid in blast_id:
            return rec
    return None


def blastn_find_region(query_fna, assembly_fna, db_dir):
    """BLAST promoter ref against assembly; return (contig, start, end, strand, pident)."""
    db = os.path.join(db_dir, "strain_db")
    run(f"makeblastdb -in {assembly_fna} -dbtype nucl -out {db} 2>/dev/null")
    out, err = run(
        f"blastn -query {query_fna} -db {db} "
        f"-outfmt '6 sseqid sstart send sstrand pident length' "
        f"-max_hsps 1 -max_target_seqs 1 "
        f"-perc_identity 60 -evalue 1e-3 -word_size 11"
    )
    if not out:
        return None
    parts = out.split("\t")
    if len(parts) < 5:
        return None
    contig = parts[0]
    s, e   = int(parts[1]), int(parts[2])
    strand = parts[3]
    pident = float(parts[4])
    return contig, min(s,e), max(s,e), strand, pident


def extract_hit(records, contig_id, start, end, strand):
    """Extract hit sequence from pre-loaded assembly records."""
    rec = match_contig(records, contig_id)
    if rec is None:
        return None
    subseq = rec.seq[start-1:end]
    if strand == "minus":
        subseq = subseq.reverse_complement()
    return str(subseq).upper()


def mafft_align(ref_seq, query_seq, tmp_dir):
    in_fa  = os.path.join(tmp_dir, "pair.fa")
    out_fa = os.path.join(tmp_dir, "aligned.fa")
    with open(in_fa, "w") as fh:
        fh.write(f">ref\n{ref_seq}\n>query\n{query_seq}\n")
    run(f"mafft --quiet --auto {in_fa} > {out_fa}")
    seqs, cur = {}, None
    for line in Path(out_fa).read_text().split("\n"):
        if line.startswith(">"):
            cur = line[1:].split()[0]; seqs[cur] = []
        elif cur:
            seqs[cur].append(line)
    if "ref" not in seqs or "query" not in seqs:
        return None, None
    return "".join(seqs["ref"]), "".join(seqs["query"])


def count_variants(ref_aln, qry_aln):
    snps = sum(1 for r,q in zip(ref_aln, qry_aln) if r!="-" and q!="-" and r!=q)
    indels = sum(1 for r,q in zip(ref_aln, qry_aln) if r=="-" or q=="-")
    return snps, indels


def analyse_target(name, ref_fna, assembly_fna, assembly_records, tmp_dir):
    result = {
        "target": name, "status": "no_hit", "pident": None,
        "snps": None, "indels": None,
        "known_variant": "none", "drug_impact": "", "notes": "",
    }

    ref_seq = "".join(
        l for l in Path(ref_fna).read_text().split("\n") if not l.startswith(">")
    ).upper()

    if not ref_seq.strip():
        result["notes"] = "Reference promoter FASTA is empty — re-run extract_promoter_refs.py"
        return result

    hit = blastn_find_region(ref_fna, assembly_fna, tmp_dir)
    if hit is None:
        result["notes"] = "BLAST: no hit (region absent or too divergent)"
        return result

    contig, start, end, strand, pident = hit
    query_seq = extract_hit(assembly_records, contig, start, end, strand)
    if query_seq is None:
        result["notes"] = f"Could not extract contig {contig}"
        return result

    result["status"] = "hit"
    result["pident"] = round(pident, 1)

    ref_aln, qry_aln = mafft_align(ref_seq, query_seq, tmp_dir)
    if ref_aln:
        result["snps"], result["indels"] = count_variants(ref_aln, qry_aln)

    # ── Check known variants ──────────────────────────────────────────────────
    variants_found = []
    if name == "mtrR":
        # Search directly in the extracted query sequence (not via reference
        # comparison), because FA1090 reference itself lacks the classic IR.
        ir_present = (MTRR_IR in query_seq) or (MTRR_IR_RC in query_seq)

        if ir_present:
            # Susceptible-type promoter: IR is intact
            variants_found.append("mtrR_promoter_IR_intact")
            result["notes"] += "13-bp IR present → susceptible mtrCDE promoter; "
        else:
            # IR absent — could be deletion (resistant) or divergent sequence
            if result["indels"] and result["indels"] > 0:
                variants_found.append("mtrR_promoter_indel")
                result["drug_impact"] = "PENICILLIN/TETRACYCLINE/MACROLIDE"
                result["notes"] += (
                    f"IR absent + {result['indels']} indel(s) vs FA1090 → "
                    "putative mtrCDE promoter structural variant; "
                )
            else:
                variants_found.append("mtrR_promoter_IR_absent")
                result["notes"] += (
                    "IR absent, no indels vs FA1090 → "
                    "likely same IR-deleted state as FA1090 reference; "
                )

    result["known_variant"] = "|".join(variants_found) if variants_found else "none"
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--assembly",  required=True)
    p.add_argument("--promoters", required=True)
    p.add_argument("--out",       required=True)
    args = p.parse_args()

    strain   = Path(args.assembly).stem
    ref_files = sorted(Path(args.promoters).glob("*_promoter_ref.fna"))
    if not ref_files:
        print(f"ERROR: no *_promoter_ref.fna in {args.promoters}")
        return

    assembly_records = load_fasta_index(args.assembly)

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for ref_fna in ref_files:
            name = ref_fna.name.replace("_promoter_ref.fna", "")
            res  = analyse_target(name, str(ref_fna), args.assembly,
                                  assembly_records, tmp)
            res["strain"] = strain
            rows.append(res)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    cols = ["strain","target","status","pident","snps","indels",
            "known_variant","drug_impact","notes"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c,"")) for c in cols) + "\n")

    print(f"  {strain}: saved → {args.out}")


if __name__ == "__main__":
    main()
