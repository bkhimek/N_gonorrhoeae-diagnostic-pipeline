#!/usr/bin/env python3
import argparse, csv
from pathlib import Path

def sniff_delim(path):
    with open(path, 'r', newline='') as f:
        sample = f.read(4096)
        try:
            return csv.Sniffer().sniff(sample, delimiters=[',','\t']).delimiter
        except Exception:
            return '\t' if '\t' in sample else ','

def main():
    ap = argparse.ArgumentParser(
        description="Select core genes with average identity >= threshold."
    )
    ap.add_argument("--input", required=True, help="Input CSV/TSV with columns incl. Gene_File (or Gene_Name) and Average_Identity(%)")
    ap.add_argument("--output", required=True, help="Output TXT with one selected gene per line")
    ap.add_argument("--thr", type=float, default=98.0, help="Identity cutoff (in percent). Default: 98.0")
    args = ap.parse_args()

    inp = Path(args.input)
    outp = Path(args.output)
    outp.parent.mkdir(parents=True, exist_ok=True)

    delim = sniff_delim(str(inp))
    selected = []

    with open(inp, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter=delim)
        for row in reader:
            # Find identity column
            ident_str = row.get("Average_Identity(%)") or row.get("Average_Identity") or ""
            try:
                identity = float(ident_str)
            except Exception:
                continue
            if identity >= args.thr:
                gene = row.get("Gene_File") or row.get("Gene_Name") or ""
                if gene:
                    selected.append(gene)

    with open(outp, "w", newline='') as outfile:
        for gene in selected:
            outfile.write(gene + "\n")

    print(f"✅ Selected {len(selected)} core genes with ≥{args.thr}% identity.")
    print(f"✅ Saved to: {outp}")

if __name__ == "__main__":
    main()
