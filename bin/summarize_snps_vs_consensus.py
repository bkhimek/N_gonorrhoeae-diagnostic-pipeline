#!/usr/bin/env python3
"""
Summarize SNPs per gene by comparing each aligned sequence to its gene consensus.
- Inputs:
  --split_dir: per-gene aligned FASTAs (same folder as above)
  --consensus_fasta: consensus per gene from build_consensus_from_split.py
- Outputs:
  --out_tsv: table with (Gene, StrainID, SNPs, Length, SNP_rate)
  --pergene_tsv: table with per-gene summary (Gene, n_strains, mean_SNPs, max_SNPs)
"""
import argparse
from pathlib import Path
from statistics import mean
from Bio import SeqIO

def snps_vs_cons(cons, seq):
    return sum(1 for a, b in zip(cons, seq) if a != b)  # gaps count as mismatches

def main():
    ap = argparse.ArgumentParser(description="Summarize SNPs per gene vs consensus.")
    ap.add_argument("--split_dir", required=True, help="Folder with per-gene aligned FASTAs (*.fasta)")
    ap.add_argument("--consensus_fasta", required=True, help="Multi-FASTA of per-gene consensus")
    ap.add_argument("--out_tsv", required=True, help="Per-strain SNP table TSV")
    ap.add_argument("--pergene_tsv", required=True, help="Per-gene summary TSV")
    args = ap.parse_args()

    split_dir = Path(args.split_dir)
    cons_map = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.consensus_fasta, "fasta")}

    rows = []
    per_gene = {}

    for fp in sorted(p for p in split_dir.iterdir() if p.suffix.lower() in (".fasta", ".fa")):
        gene = fp.stem
        cons = cons_map.get(gene)
        if not cons:
            continue
        recs = list(SeqIO.parse(str(fp), "fasta"))
        L = len(cons)
        snp_counts = []
        for r in recs:
            snp = snps_vs_cons(cons, str(r.seq))
            rows.append((gene, r.id, snp, L, (snp / L) if L else 0.0))
            snp_counts.append(snp)
        if snp_counts:
            per_gene[gene] = (len(snp_counts), mean(snp_counts), max(snp_counts))

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w") as f:
        f.write("Gene\tStrainID\tSNPs\tLength\tSNP_rate\n")
        for gene, sid, snp, L, rate in rows:
            f.write(f"{gene}\t{sid}\t{snp}\t{L}\t{rate:.6f}\n")

    with open(args.pergene_tsv, "w") as f:
        f.write("Gene\tn_strains\tmean_SNPs\tmax_SNPs\n")
        for gene in sorted(per_gene):
            n, m, mx = per_gene[gene]
            f.write(f"{gene}\t{n}\t{m:.3f}\t{mx}\n")

    print(f"[OK] Wrote SNP per-strain -> {args.out_tsv}")
    print(f"[OK] Wrote SNP per-gene  -> {args.pergene_tsv}")

if __name__ == "__main__":
    main()
