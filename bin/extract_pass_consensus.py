#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--consensus", required=True, help="Consensus FASTA (all genes)")
    ap.add_argument("--exclusivity_tsv", help="Exclusivity summary TSV (PASS/REJECT). If not given, use latest run.")
    ap.add_argument("--outdir", required=True, help="Output folder for per-gene FASTAs (PASS only)")
    args = ap.parse_args()

    # If no exclusivity TSV provided, use latest run
    excl_path = args.exclusivity_tsv
    if excl_path is None:
        runs = sorted(Path("results/runs").glob("run-*/exclusivity_summary.tsv"))
        if not runs:
            raise FileNotFoundError("No exclusivity_summary.tsv found in results/runs/*/")
        excl_path = runs[-1]  # most recent run
        print(f"[INFO] No --exclusivity_tsv provided; using latest: {excl_path}")

    excl = pd.read_csv(excl_path, sep="\t")
    pass_genes = set(excl.loc[excl["decision"]=="PASS","Gene"])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n_written = 0
    for rec in SeqIO.parse(args.consensus, "fasta"):
        if rec.id in pass_genes:
            SeqIO.write(rec, outdir / f"{rec.id}.fasta", "fasta")
            n_written += 1

    print(f"[OK] Extracted {n_written} PASS consensus FASTAs -> {outdir}")

if __name__ == "__main__":
    main()
