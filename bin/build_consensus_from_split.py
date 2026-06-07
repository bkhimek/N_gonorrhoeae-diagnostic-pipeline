#!/usr/bin/env python3
"""
Build majority-rule consensus sequences per gene from a folder of MAFFT-aligned FASTAs.
- Input: directory with *.fasta (each is a per-gene alignment)
- Optional: list of gene IDs to include (one per line)
- Output: FASTA of per-gene consensus sequences
"""
import argparse
from pathlib import Path
from collections import Counter
from Bio import SeqIO

def consensus_of_alignment(records, gap_char='-', include_gaps=False, tie_char='N'):
    seqs = [str(r.seq) for r in records]
    if not seqs:
        return ""
    L = len(seqs[0])
    cols = zip(*seqs)
    cons = []
    for col in cols:
        col = list(col)
        if not include_gaps:
            col = [c for c in col if c != gap_char]
        if not col:
            cons.append(gap_char if include_gaps else tie_char)
            continue
        counts = Counter(col)
        top, n = counts.most_common(1)[0]
        # Tie handling: if tie for top, use tie_char
        if list(counts.values()).count(n) > 1:
            cons.append(tie_char)
        else:
            cons.append(top)
    return "".join(cons)

def main():
    ap = argparse.ArgumentParser(description="Build consensus per gene from split alignments.")
    ap.add_argument("--split_dir", required=True, help="Folder with per-gene aligned FASTAs (*.fasta)")
    ap.add_argument("--out_fasta", required=True, help="Output multi-FASTA of consensus sequences")
    ap.add_argument("--genes_txt", help="Optional: file with gene names to include (one per line)")
    ap.add_argument("--include_gaps", action="store_true", help="Keep gaps in consensus")
    args = ap.parse_args()

    split_dir = Path(args.split_dir)
    out_fa = Path(args.out_fasta)
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    allow = None
    if args.genes_txt:
        with open(args.genes_txt) as f:
            allow = {ln.strip() for ln in f if ln.strip()}

    files = sorted([p for p in split_dir.iterdir() if p.suffix.lower() in (".fasta", ".fa")])

    n_written = 0
    with open(out_fa, "w") as out:
        for fp in files:
            gene = fp.stem
            if allow and gene not in allow:
                continue
            recs = list(SeqIO.parse(str(fp), "fasta"))
            if not recs:
                continue
            cons = consensus_of_alignment(recs, include_gaps=args.include_gaps)
            if not cons:
                continue
            out.write(f">{gene}\n{cons}\n")
            n_written += 1

    print(f"[OK] Wrote {n_written} consensus sequences -> {out_fa}")

if __name__ == "__main__":
    main()
