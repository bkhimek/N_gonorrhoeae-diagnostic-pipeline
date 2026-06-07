#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import subprocess

def log(msg):
    print(msg, file=sys.stdout, flush=True)

def read_core_gene_list(gpa_path: Path, min_presence: float) -> list[str]:
    """
    Decide which clusters are 'core' by presence fraction across sample columns.
    GPA columns: ['Gene','Non-unique Gene name','Annotation', <one column per sample>]
    Presence = non-empty / non-null cell in a sample column.
    """
    gpa = pd.read_csv(gpa_path)
    # sample columns are all columns except these 3:
    ann_cols = ["Gene", "Non-unique Gene name", "Annotation"]
    sample_cols = [c for c in gpa.columns if c not in ann_cols]

    if not sample_cols:
        raise RuntimeError("No sample columns found in gene_presence_absence.csv")

    # presence mask: notna AND non-empty after string conversion
    present_mask = gpa[sample_cols].notna() & gpa[sample_cols].astype(str).applymap(lambda s: len(s.strip()) > 0)
    presence_frac = present_mask.sum(axis=1) / float(len(sample_cols))

    core = gpa.loc[presence_frac >= min_presence, "Gene"].astype(str).tolist()
    return core

def load_gene_data(gene_data_path: Path) -> pd.DataFrame:
    """
    Panaroo gene_data.csv has columns like:
    ['gff_file','scaffold_name','clustering_id','annotation_id','prot_sequence',
     'dna_sequence','gene_name','description']
    We need mapping clustering_id -> gene_name (which matches GPA 'Gene').
    """
    df = pd.read_csv(gene_data_path)
    needed = {"clustering_id", "gene_name"}
    if not needed.issubset(df.columns):
        raise RuntimeError(f"{gene_data_path} missing required columns {needed}")
    # coerce to str to avoid numeric surprises
    df["clustering_id"] = df["clustering_id"].astype(str)
    df["gene_name"] = df["gene_name"].astype(str)
    return df[["clustering_id", "gene_name"]]

def load_cds_fasta(cds_path: Path) -> dict[str, str]:
    """
    combined_DNA_CDS.fasta headers are 'clustering_id' (e.g., 0_0_0).
    Return dict id -> sequence (as string).
    """
    seqs = {}
    for rec in SeqIO.parse(str(cds_path), "fasta"):
        seqs[str(rec.id)] = str(rec.seq)
    if not seqs:
        raise RuntimeError(f"No sequences read from {cds_path}")
    return seqs

def write_fasta(records: list[tuple[str, str]], path: Path):
    with path.open("w") as out:
        for rid, seq in records:
            out.write(f">{rid}\n")
            # wrap to 60 chars per line (simple wrapping)
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

def align_with_mafft(raw_fa: Path, aligned_fa: Path, threads: int):
    # mafft --thread N --auto raw.fa > aligned.fa
    cmd = ["mafft", "--thread", str(threads), "--auto", str(raw_fa)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT failed for {raw_fa}:\n{proc.stderr}")
    aligned_fa.write_text(proc.stdout)

def main():
    ap = argparse.ArgumentParser(
        description="Split Panaroo core genes into per-gene FASTAs (MAFFT-aligned). "
                    "Uses gene_presence_absence.csv + gene_data.csv + combined_DNA_CDS.fasta."
    )
    ap.add_argument("--gpa",        required=True, help="gene_presence_absence.csv")
    ap.add_argument("--gene_data",  required=True, help="gene_data.csv (from Panaroo output)")
    ap.add_argument("--cds",        required=True, help="combined_DNA_CDS.fasta")
    ap.add_argument("--outdir",     required=True, help="Output directory for per-gene FASTAs")
    ap.add_argument("--min_presence", type=float, default=1.00, help="Core threshold (fraction), e.g. 1.00 for 100%")
    ap.add_argument("--threads",      type=int,   default=8,    help="Threads for MAFFT")
    args = ap.parse_args()

    gpa_path = Path(args.gpa)
    gene_data_path = Path(args.gene_data)
    cds_path = Path(args.cds)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Choose core clusters by presence fraction
    core_genes = read_core_gene_list(gpa_path, args.min_presence)
    log(f"[INFO] core≥{args.min_presence:.2f} clusters: {len(core_genes)}")

    # 2) Load mapping clustering_id -> gene_name
    gd = load_gene_data(gene_data_path)

    # 3) Load CDS sequences by clustering_id
    cds = load_cds_fasta(cds_path)

    wrote = 0
    skipped = 0
    for gene in core_genes:
        # all ids that belong to this cluster name
        ids = gd.loc[gd["gene_name"] == gene, "clustering_id"].astype(str).tolist()
        # map to sequences (some ids may be absent if filtered earlier)
        recs = [(cid, cds[cid]) for cid in ids if cid in cds]

        raw_path = outdir / f"{gene}.raw.fa"
        aln_path = outdir / f"{gene}.fasta"

        if len(recs) < 2:
            # write raw (for debugging) but skip MAFFT
            if recs:
                write_fasta(recs, raw_path)
            skipped += 1
            continue

        write_fasta(recs, raw_path)
        align_with_mafft(raw_path, aln_path, args.threads)
        wrote += 1

    log(f"[OK] wrote {wrote} aligned per-gene FASTAs to: {outdir}")
    log(f"[NOTE] skipped {skipped} clusters with <2 mapped sequences")

if __name__ == "__main__":
    main()
