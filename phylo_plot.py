#!/usr/bin/env python3
"""
phylo_plot.py — Annotated N. gonorrhoeae phylogenetic tree (interactive HTML)

Reads:
  - parsnp.tree   (Newick from Parsnp)
  - AMR TSV dir   (AMRFinderPlus per-strain results)

Outputs:
  - phylo_amr_tree.html  (interactive, opens in any browser)

Usage:
    python phylo_plot.py \
        --tree ~/n_gonorrhoeae_project/phylogeny/parsnp.tree \
        --amr  ~/n_gonorrhoeae_project/results/amr/ \
        --out  ~/n_gonorrhoeae_project/phylogeny/
"""

import argparse
import glob
import os

import pandas as pd
import plotly.graph_objects as go
from Bio import Phylo

# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True, help="parsnp.tree (Newick)")
    p.add_argument("--amr",  required=True, help="Directory of AMRFinder TSVs")
    p.add_argument("--out",  default=".",   help="Output directory")
    return p.parse_args()

# ─── AMR burden ───────────────────────────────────────────────────────────────

KNOWN_SUFFIXES = ("_amrfinder", "_amr", "_amrfinderplus", "_amr_report")

def load_amr_burden(amr_dir):
    files = glob.glob(os.path.join(amr_dir, "*.tsv"))
    burden, top_class = {}, {}
    for f in files:
        strain = os.path.basename(f).rsplit(".", 1)[0]
        for sfx in KNOWN_SUFFIXES:
            strain = strain.replace(sfx, "")
        try:
            df = pd.read_csv(f, sep="\t", low_memory=False)
            df.columns = df.columns.str.strip()
            if "Element type" in df.columns:
                df = df[df["Element type"].isin(["AMR", "POINTP"])]
            burden[strain] = df["Element symbol"].nunique() if "Element symbol" in df.columns else 0
            if "Class" in df.columns and not df.empty:
                top_class[strain] = df["Class"].value_counts().idxmax()
            else:
                top_class[strain] = "Unknown"
        except Exception:
            burden[strain] = 0
            top_class[strain] = "Unknown"
    return burden, top_class

# ─── Tree layout ──────────────────────────────────────────────────────────────

def compute_layout(tree):
    tips = tree.get_terminals()
    tip_order = {t.name: i for i, t in enumerate(tips)}
    node_x, node_y = {}, {}

    def assign_x(clade, px):
        bl = clade.branch_length or 0.0
        node_x[id(clade)] = px + bl
        for c in clade.clades:
            assign_x(c, px + bl)

    tree.root.branch_length = 0.0
    assign_x(tree.root, 0.0)

    def assign_y(clade):
        if clade.is_terminal():
            y = tip_order[clade.name]
        else:
            child_ys = [assign_y(c) for c in clade.clades]
            y = (min(child_ys) + max(child_ys)) / 2.0
        node_y[id(clade)] = y
        return y

    assign_y(tree.root)
    return node_x, node_y, tips


def get_line_coords(tree, node_x, node_y):
    lx, ly = [], []

    def dfs(clade):
        x, y = node_x[id(clade)], node_y[id(clade)]
        if not clade.is_terminal():
            child_ys = [node_y[id(c)] for c in clade.clades]
            lx.extend([x, x, None])
            ly.extend([min(child_ys), max(child_ys), None])
            for c in clade.clades:
                lx.extend([x, node_x[id(c)], None])
                ly.extend([node_y[id(c)], node_y[id(c)], None])
                dfs(c)

    dfs(tree.root)
    return lx, ly

# ─── Build and save figure ────────────────────────────────────────────────────

def make_figure(tree_file, amr_dir, out_dir):
    print("Parsing tree ...")
    tree = Phylo.read(tree_file, "newick")

    print("Loading AMR burden ...")
    burden, top_class = load_amr_burden(amr_dir)

    print("Computing layout ...")
    node_x, node_y, tips = compute_layout(tree)
    lx, ly = get_line_coords(tree, node_x, node_y)

    n = len(tips)
    max_b = max(burden.values()) if burden else 1

    # Branch lines
    line_trace = go.Scatter(
        x=lx, y=ly, mode="lines",
        line=dict(color="#aaaaaa", width=0.5),
        hoverinfo="skip", showlegend=False,
    )

    # Tip markers
    tip_x, tip_y, tip_text, tip_color = [], [], [], []
    for t in tips:
        strain = t.name.replace(".fna", "")
        b = burden.get(strain, 0)
        tc = top_class.get(strain, "Unknown")
        tip_x.append(node_x[id(t)])
        tip_y.append(node_y[id(t)])
        tip_text.append(f"<b>{strain}</b><br>AMR elements: {b}<br>Top class: {tc}")
        tip_color.append(b)

    tip_trace = go.Scatter(
        x=tip_x, y=tip_y, mode="markers",
        marker=dict(
            color=tip_color, colorscale="YlOrRd", size=6,
            cmin=0, cmax=max_b,
            colorbar=dict(
                title=dict(text="AMR gene<br>burden", font=dict(size=11)),
                thickness=14, len=0.5, yanchor="middle", y=0.5,
            ),
            line=dict(width=0.3, color="#666666"),
        ),
        text=tip_text,
        hovertemplate="%{text}<extra></extra>",
        showlegend=False,
    )

    fig = go.Figure([line_trace, tip_trace])
    fig.update_layout(
        title=dict(
            text=(f"<b>N. gonorrhoeae</b> core-genome phylogeny "
                  f"annotated with AMR burden ({n} strains)<br>"
                  "<sup>Parsnp v2.1.5 · FA1909 reference · AMRFinderPlus v4.0.23</sup>"),
            font=dict(size=14),
        ),
        xaxis=dict(title="SNP distance from root", showgrid=False),
        yaxis=dict(visible=False),
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=max(700, n * 5),
        margin=dict(l=20, r=100, t=80, b=50),
        hovermode="closest",
    )

    out_path = os.path.join(out_dir, "phylo_amr_tree.html")
    fig.write_html(out_path, include_plotlyjs="cdn")
    print(f"Saved → {out_path}")
    print(f"Open in browser: file://{os.path.abspath(out_path)}")


if __name__ == "__main__":
    args = parse_args()
    make_figure(args.tree, args.amr, args.out)
