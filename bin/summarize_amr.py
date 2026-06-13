#!/usr/bin/env python3
"""
summarize_amr.py
Merges per-gene mutation TSVs + AMRFinderPlus TSVs into a wide
resistance_profiles.tsv and an HTML summary report.
"""
import argparse, csv, sys
from pathlib import Path
from collections import defaultdict

def load_mutations(mutations_dir):
    """Returns dict: strain_id -> {label -> is_resistance bool}"""
    data = defaultdict(dict)
    for tsv in Path(mutations_dir).glob("*_mutations.tsv"):
        with open(tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid   = row["strain_id"]
                label = row["label"]
                data[sid][label] = row["is_resistance"].lower() == "true"
    return data

def load_amrfinder(amrfinder_dir):
    """Returns dict: strain_id -> list of element_symbol strings"""
    data = defaultdict(list)
    for tsv in Path(amrfinder_dir).glob("*_amrfinder.tsv"):
        with open(tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = row.get("Name", "")
                sym = row.get("Element symbol", "")
                if sid and sym:
                    data[sid].append(sym)
    return data

PROFILE_COLS = [
    "strain_id",
    "GyrA_S91", "GyrA_D95",
    "ParC_D86", "ParC_S87", "ParC_E91",
    "PenA_mosaic", "PenA_A345G", "PenA_D346ins", "PenA_A501V", "PenA_G543",
    "MtrR_A39T",
    "PorB_G120", "PorB_A121",
    "amrfinder_hits",
    "predicted_cipro_r",
    "predicted_ceph_r",
]

def predict_cipro(row):
    """Fluoroquinolone resistance: GyrA S91 OR D95 mutation."""
    return row.get("GyrA_S91", False) or row.get("GyrA_D95", False)

def predict_ceph(row):
    """Reduced cephalosporin susceptibility: penA mosaic allele."""
    return row.get("PenA_mosaic", False)

def build_profiles(mut_data, amr_data):
    all_strains = sorted(set(list(mut_data.keys()) + list(amr_data.keys())))
    rows = []
    for sid in all_strains:
        m = mut_data.get(sid, {})
        row = {"strain_id": sid}
        for col in PROFILE_COLS[1:-3]:   # mutation columns
            row[col] = m.get(col, False)
        row["amrfinder_hits"] = "|".join(amr_data.get(sid, [])) or "none"
        row["predicted_cipro_r"] = predict_cipro(row)
        row["predicted_ceph_r"]  = predict_ceph(row)
        rows.append(row)
    return rows

def write_tsv(rows, outpath):
    with open(outpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PROFILE_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def write_html(rows, outpath):
    cipro_r = sum(1 for r in rows if r["predicted_cipro_r"])
    ceph_r  = sum(1 for r in rows if r["predicted_ceph_r"])
    total   = len(rows)
    header  = "\n".join(f"<th>{c}</th>" for c in PROFILE_COLS)
    trows   = ""
    for r in rows:
        cells = ""
        for c in PROFILE_COLS:
            val = r[c]
            style = ""
            if val is True:
                style = ' style="background:#ffcccc"'
            elif val is False and c not in ("strain_id","amrfinder_hits",
                                            "predicted_cipro_r","predicted_ceph_r"):
                style = ' style="background:#ccffcc"'
            cells += f"<td{style}>{val}</td>"
        trows += f"<tr>{cells}</tr>\n"
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>N. gonorrhoeae AMR Summary</title>
<style>body{{font-family:monospace;font-size:12px}}
table{{border-collapse:collapse}}
th,td{{border:1px solid #ccc;padding:3px 6px;white-space:nowrap}}
th{{background:#333;color:#fff}}</style></head><body>
<h2>N. gonorrhoeae AMR Resistance Profile</h2>
<p><b>Total strains:</b> {total} &nbsp;|&nbsp;
<b>Predicted ciprofloxacin R:</b> {cipro_r} ({100*cipro_r//total if total else 0}%) &nbsp;|&nbsp;
<b>Predicted reduced ceph susceptibility:</b> {ceph_r} ({100*ceph_r//total if total else 0}%)</p>
<table><thead><tr>{header}</tr></thead><tbody>
{trows}</tbody></table></body></html>"""
    with open(outpath, "w") as fh:
        fh.write(html)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mutations-dir",  required=True)
    p.add_argument("--amrfinder-dir",  required=True)
    p.add_argument("--output-tsv",     required=True)
    p.add_argument("--output-html",    required=True)
    args = p.parse_args()

    mut_data = load_mutations(args.mutations_dir)
    amr_data = load_amrfinder(args.amrfinder_dir)
    rows     = build_profiles(mut_data, amr_data)

    write_tsv(rows,  args.output_tsv)
    write_html(rows, args.output_html)

    print(f"INFO: {len(rows)} strains → {args.output_tsv}, {args.output_html}",
          file=sys.stderr)

if __name__ == "__main__":
    main()
