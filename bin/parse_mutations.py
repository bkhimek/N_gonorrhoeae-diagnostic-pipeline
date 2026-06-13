#!/usr/bin/env python3
import argparse, sys
from pathlib import Path

KNOWN_MUTATIONS = {
    "gyrA": {
        91: {"ref": "S", "resistance": {"F","Y"},            "label": "GyrA_S91"},
        95: {"ref": "D", "resistance": {"G","N","A","H"},    "label": "GyrA_D95"},
    },
    "parC": {
        86: {"ref": "D", "resistance": {"N","Y"},            "label": "ParC_D86"},
        87: {"ref": "S", "resistance": {"R","N","I","T"},    "label": "ParC_S87"},
        91: {"ref": "E", "resistance": {"K","G"},            "label": "ParC_E91"},
    },
    "penA": {
        345: {"ref": "A", "resistance": {"G"},               "label": "PenA_A345G"},
        346: {"ref": "D", "resistance": {"DD"},              "label": "PenA_D346ins"},
        501: {"ref": "A", "resistance": {"V"},               "label": "PenA_A501V"},
        543: {"ref": "G", "resistance": {"S","D"},           "label": "PenA_G543"},
    },
    "mtrR": {
        39:  {"ref": "A", "resistance": {"T"},               "label": "MtrR_A39T"},
    },
    "porB": {
        120: {"ref": "G", "resistance": {"D","K","N"},       "label": "PorB_G120"},
        121: {"ref": "A", "resistance": {"D","P","N"},       "label": "PorB_A121"},
    },
}
PENA_MOSAIC_THRESHOLD = 5

def parse_fasta(path):
    seqs, header, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(chunks).upper()))
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        seqs.append((header, "".join(chunks).upper()))
    return seqs

def build_ref_col_map(ref_seq):
    ref_col, ref_pos = {}, 0
    for col, aa in enumerate(ref_seq):
        if aa != "-":
            ref_pos += 1
            ref_col[ref_pos] = col
    return ref_col

def call_known_positions(sid, gene, ref_seq, qry_seq, ref_col_map, known):
    rows = []
    for ref_pos, meta in known.items():
        col = ref_col_map.get(ref_pos)
        if col is None: continue
        ref_aa, qry_aa = ref_seq[col], qry_seq[col]
        if qry_aa == "-": continue
        rows.append({"strain_id": sid, "gene": gene, "label": meta["label"],
                     "ref_position": ref_pos, "ref_aa": ref_aa, "alt_aa": qry_aa,
                     "is_known_site": True,
                     "is_resistance": qry_aa in meta["resistance"], "note": ""})
    return rows

def call_qrdr_scan(sid, gene, ref_seq, qry_seq, ref_col_map, start, end, known_pos):
    rows = []
    for ref_pos in range(start, end + 1):
        if ref_pos in known_pos: continue
        col = ref_col_map.get(ref_pos)
        if col is None: continue
        ref_aa, qry_aa = ref_seq[col], qry_seq[col]
        if qry_aa == "-" or qry_aa == ref_aa: continue
        rows.append({"strain_id": sid, "gene": gene,
                     "label": f"{gene}_{ref_aa}{ref_pos}",
                     "ref_position": ref_pos, "ref_aa": ref_aa, "alt_aa": qry_aa,
                     "is_known_site": False, "is_resistance": False,
                     "note": "novel_qrdr_variant"})
    return rows

def call_pena_mosaic(sid, ref_seq, qry_seq):
    n = sum(1 for r,q in zip(ref_seq,qry_seq) if r!="-" and q!="-" and r!=q)
    return {"strain_id": sid, "gene": "penA", "label": "PenA_mosaic",
            "ref_position": -1, "ref_aa": "", "alt_aa": "",
            "is_known_site": True, "is_resistance": n >= PENA_MOSAIC_THRESHOLD,
            "note": f"aa_substitutions={n}"}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--alignment",  required=True)
    p.add_argument("--gene",       required=True)
    p.add_argument("--reference",  required=True)
    p.add_argument("--qrdr-start", type=int)
    p.add_argument("--qrdr-end",   type=int)
    p.add_argument("--output",     required=True)
    args = p.parse_args()

    gene = args.gene
    aln  = parse_fasta(Path(args.alignment))
    if not aln:
        sys.exit(f"ERROR: empty alignment: {args.alignment}")

    refs = [(h,s) for h,s in aln if args.reference in h]
    if not refs:
        sys.exit(f"ERROR: reference '{args.reference}' not found")
    _, ref_seq = refs[0]

    ref_col_map  = build_ref_col_map(ref_seq)
    known        = KNOWN_MUTATIONS.get(gene, {})
    do_qrdr      = args.qrdr_start is not None and args.qrdr_end is not None
    COLS = ["strain_id","gene","label","ref_position",
            "ref_aa","alt_aa","is_known_site","is_resistance","note"]
    rows = []

    for header, seq in aln:
        if args.reference in header: continue
        sid = header.split("__")[0]
        rows.extend(call_known_positions(sid, gene, ref_seq, seq, ref_col_map, known))
        if do_qrdr:
            rows.extend(call_qrdr_scan(sid, gene, ref_seq, seq, ref_col_map,
                                       args.qrdr_start, args.qrdr_end, set(known)))
        if gene == "penA":
            rows.append(call_pena_mosaic(sid, ref_seq, seq))

    with open(args.output, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in COLS) + "\n")

    print(f"INFO: {gene} — {len(rows)} rows, {len(aln)-1} strains → {args.output}",
          file=sys.stderr)

if __name__ == "__main__":
    main()
