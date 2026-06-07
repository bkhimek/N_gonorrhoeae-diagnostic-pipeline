#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from pathlib import Path

def clean_neighbor(s: str) -> str:
    """Extract species name from BLAST stitle field."""
    if not isinstance(s, str):
        return "unknown"
    # Drop accession prefix like 'NZ_LR594047.1 '
    s = re.sub(r"^[A-Z0-9_.]+\s+", "", s).strip()
    # Cut off trailing descriptors
    s = re.sub(r"\b(strain|chromosome|complete genome|contig|isolate|genome)\b.*",
               "", s, flags=re.I).strip()
    s = re.sub(r"\s+", " ", s)
    return s if s else "unknown"

def main():
    ap = argparse.ArgumentParser(description="Summarize BLAST hits for exclusivity check")
    ap.add_argument("--blast_tsv", required=True, help="BLAST outfmt 6 (qseqid..qcovs+sscinames,stitle)")
    ap.add_argument("--out_tsv",   required=True, help="Output summary TSV")
    ap.add_argument("--min_pident", type=float, default=85.0, help="Min %% identity to call REJECT")
    ap.add_argument("--min_qcovs",  type=float, default=80.0, help="Min %% query coverage to call REJECT")
    args = ap.parse_args()

    cols = [
        "qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore",
        "qlen","slen","qcovs","sscinames","staxids","stitle"
    ]

    p = Path(args.blast_tsv)
    if not p.exists() or p.stat().st_size == 0:
        pd.DataFrame(columns=["Gene","neighbor","best_pident","best_qcovs","decision"])\
          .to_csv(args.out_tsv, sep="\t", index=False)
        print(f"[OK] No hits found; all genes PASS -> {args.out_tsv}")
        return

    # Important: header=None so first line is data, not a header row
    df = pd.read_csv(args.blast_tsv, sep="\t", header=None, names=cols)

    if df.empty:
        pd.DataFrame(columns=["Gene","neighbor","best_pident","best_qcovs","decision"])\
          .to_csv(args.out_tsv, sep="\t", index=False)
        print(f"[OK] No hits found; all genes PASS -> {args.out_tsv}")
        return

    # Clean neighbor name from stitle (sscinames is often NaN here)
    df["neighbor"] = df["stitle"].map(clean_neighbor)

    # Pick best hit per gene: highest identity, then coverage, then bitscore
    best = (df.sort_values(["qseqid","pident","qcovs","bitscore"],
                           ascending=[True, False, False, False])
              .groupby("qseqid", as_index=False)
              .first())

    # Prepare output table
    out = best[["qseqid","neighbor","pident","qcovs"]].copy()
    out.rename(columns={"qseqid":"Gene", "pident":"best_pident", "qcovs":"best_qcovs"}, inplace=True)

    # Ensure numeric (guard against accidental strings)
    out["best_pident"] = pd.to_numeric(out["best_pident"], errors="coerce")
    out["best_qcovs"]  = pd.to_numeric(out["best_qcovs"],  errors="coerce")

    def decide(row):
        return "REJECT" if (row["best_pident"] >= args.min_pident and row["best_qcovs"] >= args.min_qcovs) else "PASS"

    out["decision"] = out.apply(decide, axis=1)

    out.sort_values("best_pident", inplace=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[OK] Exclusivity summary -> {args.out_tsv}")
    print(out["decision"].value_counts(dropna=False).to_string())

if __name__ == "__main__":
    main()

