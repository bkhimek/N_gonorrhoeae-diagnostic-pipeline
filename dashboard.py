#!/usr/bin/env python3
"""
AMR Dashboard — N. gonorrhoeae diagnostic pipeline
Visualises AMRFinderPlus results + core-genome phylogeny.

Usage:
    python dashboard.py --data ~/n_gonorrhoeae_project/results/amr/
    python dashboard.py --data ~/n_gonorrhoeae_project/results/amr/ \
                        --tree ~/n_gonorrhoeae_project/phylogeny/parsnp.tree
"""

import argparse
import glob
import os
import sys

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, Input, Output, dcc, html

# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data", default="results/amr",
                   help="Directory containing AMRFinderPlus .tsv files")
    p.add_argument("--tree", default=None,
                   help="parsnp.tree Newick file (optional, enables phylogeny tab)")
    p.add_argument("--port", type=int, default=8050)
    return p.parse_args()

# ─── AMR data loading ─────────────────────────────────────────────────────────

KNOWN_SUFFIXES = ("_amrfinder", "_amr", "_amrfinderplus", "_amr_report")

def load_results(data_dir):
    files = glob.glob(os.path.join(data_dir, "*.tsv"))
    if not files:
        sys.exit(f"ERROR: no .tsv files found in {data_dir!r}")
    frames = []
    for f in sorted(files):
        strain = os.path.basename(f).rsplit(".", 1)[0]
        for sfx in KNOWN_SUFFIXES:
            strain = strain.replace(sfx, "")
        try:
            df = pd.read_csv(f, sep="\t", low_memory=False)
            df["Strain"] = strain
            frames.append(df)
        except Exception as exc:
            print(f"  Warning: skipping {os.path.basename(f)}: {exc}")
    combined = pd.concat(frames, ignore_index=True)
    combined.columns = combined.columns.str.strip()
    if "Element type" in combined.columns:
        combined = combined[combined["Element type"].isin(["AMR", "POINTP"])]
    return combined


def build_tables(df, n_strains):
    gene_prev = (
        df.groupby(["Element symbol", "Class"])["Strain"]
        .nunique().reset_index().rename(columns={"Strain": "N"})
    )
    gene_prev["Pct"] = (gene_prev["N"] / n_strains * 100).round(1)
    gene_prev.sort_values("Pct", ascending=False, inplace=True)
    gene_prev.reset_index(drop=True, inplace=True)

    class_prev = (
        df.groupby("Class")["Strain"]
        .nunique().reset_index().rename(columns={"Strain": "N"})
    )
    class_prev["Pct"] = (class_prev["N"] / n_strains * 100).round(1)
    class_prev.sort_values("Pct", inplace=True)

    pa = (
        df[["Strain", "Element symbol"]].drop_duplicates()
        .assign(v=1)
        .pivot(index="Strain", columns="Element symbol", values="v")
        .fillna(0).astype(int)
    )
    pa = pa.loc[pa.sum(axis=1).sort_values(ascending=False).index]

    burden = df.groupby("Strain")["Element symbol"].nunique().to_dict()
    return gene_prev, class_prev, pa, burden

# ─── Colours ──────────────────────────────────────────────────────────────────

CLASS_COLOURS = {
    "BETA-LACTAM":                        "#e45756",
    "TETRACYCLINE":                       "#f58518",
    "AMINOGLYCOSIDE":                     "#4c78a8",
    "FLUOROQUINOLONE":                    "#72b7b2",
    "MACROLIDE":                          "#54a24b",
    "SULFONAMIDE":                        "#b279a2",
    "TRIMETHOPRIM":                       "#ff9da6",
    "PHENICOL":                           "#9ecae9",
    "RIFAMYCIN":                          "#fdae6b",
    "QUINOLONE":                          "#e7ba52",
    "BETA-LACTAM/QUINOLONE/TETRACYCLINE": "#d4a0a0",
    "BETA-LACTAM/MACROLIDE/TETRACYCLINE": "#7b6ea6",
}

# ─── Figures: AMR tabs ────────────────────────────────────────────────────────

def fig_prevalence(gene_prev, top_n, sel_classes):
    data = gene_prev.copy()
    if sel_classes:
        data = data[data["Class"].isin(sel_classes)]
    data = data.head(int(top_n))
    if data.empty:
        return go.Figure().update_layout(title="No data for selected filters")
    fig = px.bar(data, x="Element symbol", y="Pct", color="Class",
                 color_discrete_map=CLASS_COLOURS,
                 hover_data={"N": True, "Pct": True, "Class": True},
                 labels={"Element symbol": "Resistance element", "Pct": "Prevalence (%)"},
                 title=f"Top {len(data)} resistance elements by prevalence")
    fig.update_layout(xaxis_tickangle=-40, legend_title="Drug class",
                      plot_bgcolor="white", paper_bgcolor="white",
                      yaxis=dict(gridcolor="#eeeeee", range=[0, 108]),
                      margin=dict(b=130))
    return fig


def fig_class_summary(class_prev, n_strains):
    fig = px.bar(class_prev, x="Pct", y="Class", orientation="h",
                 color="Class", color_discrete_map=CLASS_COLOURS, text="N",
                 labels={"Pct": "% of strains", "Class": "Drug class"},
                 title=f"Drug-class resistance prevalence  ({n_strains} strains)")
    fig.update_traces(texttemplate=" %{text} strains", textposition="outside")
    fig.update_layout(showlegend=False, plot_bgcolor="white", paper_bgcolor="white",
                      xaxis=dict(gridcolor="#eeeeee", range=[0, 115]), height=420)
    return fig


def fig_heatmap(pa, top_n):
    top_genes = pa.sum().sort_values(ascending=False).head(int(top_n)).index
    sub = pa[top_genes]
    row_h = max(6, min(18, 900 // max(len(sub), 1)))
    fig = px.imshow(sub, color_continuous_scale=["#f0f4f8", "#2166ac"], aspect="auto",
                    labels={"x": "Resistance element", "y": "Strain", "color": "Present"},
                    title=f"Presence/absence — top {len(top_genes)} elements across {len(sub)} strains")
    fig.update_layout(coloraxis_showscale=False, xaxis_tickangle=-50,
                      height=max(500, len(sub) * row_h), margin=dict(l=120, b=120))
    return fig

# ─── Figure: phylogeny tab ────────────────────────────────────────────────────

def build_tree_traces(tree_file, burden):
    """Parse Newick tree and return Plotly traces for a rectangular phylogram."""
    try:
        from Bio import Phylo
    except ImportError:
        return None, "BioPython not installed. Run: conda install -c conda-forge biopython"

    tree = Phylo.read(tree_file, "newick")
    tips = tree.get_terminals()
    tip_order = {t.name: i for i, t in enumerate(tips)}
    node_x, node_y = {}, {}

    # Top-down: x = cumulative branch length
    def assign_x(clade, px):
        bl = clade.branch_length or 0.0
        node_x[id(clade)] = px + bl
        for c in clade.clades:
            assign_x(c, px + bl)

    tree.root.branch_length = 0.0
    assign_x(tree.root, 0.0)

    # Bottom-up: y = tip index or midpoint of children
    def assign_y(clade):
        if clade.is_terminal():
            y = tip_order[clade.name]
        else:
            y = (min(assign_y(c) for c in clade.clades) +
                 max(node_y[id(c)] for c in clade.clades)) / 2.0
        node_y[id(clade)] = y
        return y

    assign_y(tree.root)

    # Branch lines
    lx, ly = [], []

    def collect_lines(clade):
        x, y = node_x[id(clade)], node_y[id(clade)]
        if not clade.is_terminal():
            child_ys = [node_y[id(c)] for c in clade.clades]
            lx.extend([x, x, None])
            ly.extend([min(child_ys), max(child_ys), None])
            for c in clade.clades:
                lx.extend([x, node_x[id(c)], None])
                ly.extend([node_y[id(c)], node_y[id(c)], None])
                collect_lines(c)

    collect_lines(tree.root)

    line_trace = go.Scatter(x=lx, y=ly, mode="lines",
                            line=dict(color="#999999", width=0.5),
                            hoverinfo="skip", showlegend=False)

    # Tip scatter coloured by burden
    max_b = max(burden.values()) if burden else 1
    tip_x, tip_y, tip_text, tip_color = [], [], [], []
    for t in tips:
        strain = t.name.replace(".fna", "")
        b = burden.get(strain, 0)
        tip_x.append(node_x[id(t)])
        tip_y.append(node_y[id(t)])
        tip_text.append(f"{strain}<br>AMR elements: {b}")
        tip_color.append(b)

    tip_trace = go.Scatter(
        x=tip_x, y=tip_y, mode="markers",
        marker=dict(color=tip_color, colorscale="YlOrRd", size=5,
                    cmin=0, cmax=max_b,
                    colorbar=dict(title="AMR burden", thickness=12, len=0.4,
                                  yanchor="bottom", y=0.05)),
        text=tip_text, hovertemplate="%{text}<extra></extra>",
        showlegend=False,
    )

    n = len(tips)
    fig = go.Figure([line_trace, tip_trace])
    fig.update_layout(
        title=f"Core-genome phylogeny annotated with AMR burden ({n} strains)",
        xaxis_title="SNP distance from root",
        yaxis=dict(visible=False),
        plot_bgcolor="white", paper_bgcolor="white",
        height=max(600, n * 5),
        margin=dict(l=20, r=80, t=50, b=40),
        hovermode="closest",
    )
    return fig, None

# ─── Layout ───────────────────────────────────────────────────────────────────

HEADER_STYLE = {"backgroundColor": "#2c3e50", "color": "white", "padding": "16px 28px"}
CARD = {"border": "1px solid #dde3ea", "borderRadius": "8px",
        "padding": "16px 20px", "marginBottom": "16px", "backgroundColor": "white"}
CTRL_ROW = {"display": "flex", "alignItems": "flex-end",
            "gap": "32px", "marginBottom": "12px", "flexWrap": "wrap"}


def tab_prevalence(gene_prev, all_classes):
    n, smax = len(gene_prev), min(len(gene_prev), 50)
    return html.Div([
        html.Div([
            html.Div([
                html.Label("Top N elements", style={"fontWeight": "600"}),
                dcc.Slider(10, smax, step=5, value=min(20, n), id="top-n-slider",
                           marks={i: str(i) for i in range(10, smax + 1, 10)},
                           tooltip={"placement": "bottom", "always_visible": True}),
            ], style={"flex": "0 0 320px"}),
            html.Div([
                html.Label("Filter by drug class", style={"fontWeight": "600"}),
                dcc.Dropdown(all_classes, multi=True, id="class-filter",
                             placeholder="All classes", style={"minWidth": "280px"}),
            ], style={"flex": "1"}),
        ], style=CTRL_ROW),
        dcc.Graph(id="prev-chart", config={"displayModeBar": False}),
    ], style=CARD)


def tab_class(class_prev, n_strains):
    return html.Div([
        dcc.Graph(figure=fig_class_summary(class_prev, n_strains),
                  config={"displayModeBar": False}),
    ], style=CARD)


def tab_heatmap(pa):
    n, smax = len(pa.columns), min(len(pa.columns), 50)
    return html.Div([
        html.Div([
            html.Label("Top N elements in heatmap", style={"fontWeight": "600"}),
            dcc.Slider(10, smax, step=5, value=min(30, n), id="heat-n-slider",
                       marks={i: str(i) for i in range(10, smax + 1, 10)},
                       tooltip={"placement": "bottom", "always_visible": True}),
        ], style={"maxWidth": "380px", "marginBottom": "12px"}),
        dcc.Graph(id="heat-chart", config={"displayModeBar": False}),
    ], style=CARD)


def tab_phylo(tree_fig, err_msg):
    if err_msg:
        return html.Div([
            html.P(err_msg, style={"color": "#c0392b", "padding": "20px"}),
        ], style=CARD)
    return html.Div([
        html.P("Hover over tip dots to see strain name and AMR burden. "
               "Dots coloured by number of unique resistance elements (white=0, red=max).",
               style={"fontSize": "0.85rem", "color": "#555", "marginBottom": "8px"}),
        dcc.Graph(figure=tree_fig, config={"displayModeBar": True},
                  style={"overflowY": "auto"}),
    ], style=CARD)


def make_layout(gene_prev, class_prev, pa, n_strains, all_classes, has_tree):
    tabs = [
        dcc.Tab(label="Gene prevalence",    value="tab-prev"),
        dcc.Tab(label="Drug class summary", value="tab-class"),
        dcc.Tab(label="Strain heatmap",     value="tab-heat"),
    ]
    if has_tree:
        tabs.append(dcc.Tab(label="Phylogeny", value="tab-phylo"))

    return html.Div([
        html.Div([
            html.H2("N. gonorrhoeae AMR Dashboard",
                    style={"margin": 0, "fontSize": "1.4rem"}),
            html.P(f"{n_strains} strains · {len(gene_prev)} resistance elements · "
                   f"{len(all_classes)} drug classes",
                   style={"margin": "4px 0 0", "opacity": 0.75, "fontSize": "0.88rem"}),
        ], style=HEADER_STYLE),
        html.Div([
            dcc.Tabs(id="tabs", value="tab-prev", style={"marginBottom": "12px"},
                     children=tabs),
            html.Div(id="tab-content"),
        ], style={"maxWidth": "1400px", "margin": "0 auto", "padding": "16px 20px"}),
    ], style={"fontFamily": "Arial, sans-serif", "minHeight": "100vh",
              "backgroundColor": "#f4f6f9"})

# ─── App ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    print(f"Loading AMRFinder results from {args.data!r} ...")
    df = load_results(args.data)
    n_strains = df["Strain"].nunique()
    print(f"  {len(df):,} hits across {n_strains} strains.")

    gene_prev, class_prev, pa, burden = build_tables(df, n_strains)
    all_classes = sorted(df["Class"].dropna().unique())

    # Phylogeny (optional)
    tree_fig, tree_err = None, None
    has_tree = args.tree and os.path.isfile(args.tree)
    if has_tree:
        print(f"Building phylogeny traces from {args.tree!r} ...")
        tree_fig, tree_err = build_tree_traces(args.tree, burden)
        if tree_err:
            print(f"  Warning: {tree_err}")
        else:
            print("  Phylogeny ready.")

    app = Dash(__name__, title="N. gonorrhoeae AMR")
    app.layout = make_layout(gene_prev, class_prev, pa, n_strains, all_classes, has_tree)

    @app.callback(Output("tab-content", "children"), Input("tabs", "value"))
    def render_tab(tab):
        if tab == "tab-prev":   return tab_prevalence(gene_prev, all_classes)
        if tab == "tab-class":  return tab_class(class_prev, n_strains)
        if tab == "tab-heat":   return tab_heatmap(pa)
        if tab == "tab-phylo":  return tab_phylo(tree_fig, tree_err)

    @app.callback(Output("prev-chart", "figure"),
                  Input("top-n-slider", "value"), Input("class-filter", "value"))
    def update_prev(top_n, sel_classes):
        return fig_prevalence(gene_prev, top_n or 20, sel_classes or [])

    @app.callback(Output("heat-chart", "figure"), Input("heat-n-slider", "value"))
    def update_heat(top_n):
        return fig_heatmap(pa, top_n or 30)

    print(f"\n  Dashboard running at http://localhost:{args.port}/\n")
    app.run(debug=False, port=args.port)


if __name__ == "__main__":
    main()
