import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from data import load_hgt1, load_hgt2, dataset_sidebar
from config import TAX_LEVELS, TAX_LEVEL_NAMES

st.set_page_config(page_title="Transfer Flows", layout="wide")
st.title("Transfer Flows")

data_dir, assembly = dataset_sidebar()
st.markdown(
    "Sankey diagram showing the flow of genetic material between taxonomic groups. "
    "**Method 1**: region taxonomy (donor lineage) → contig taxonomy (recipient organism). "
    "**Method 2**: contig1 taxonomy → contig2 taxonomy for shared ORF pairs."
)

# ── Controls ──────────────────────────────────────────────────────────────────
col1, col2, col3, col4 = st.columns([2, 2, 2, 2])
with col1:
    method = st.radio("Method", ["HGT1 (Taxonomy Inference)", "HGT2 (Shared ORFs)"], key="method")
with col2:
    lvl = st.selectbox(
        "Taxonomic Level",
        options=TAX_LEVELS,
        format_func=lambda x: TAX_LEVEL_NAMES[x],
        index=1,
        key="lvl",
    )
with col3:
    min_count = st.number_input("Minimum transfer count", min_value=1, value=5, step=1)
with col4:
    if method.startswith("HGT1"):
        det_levels = st.multiselect(
            "Filter by detection level",
            options=TAX_LEVELS,
            default=TAX_LEVELS,
            format_func=lambda x: TAX_LEVEL_NAMES[x],
        )
    else:
        det_levels = None

# ── Build pairs dataframe ─────────────────────────────────────────────────────
with st.spinner("Loading..."):
    if method.startswith("HGT1"):
        df = load_hgt1(data_dir, assembly)
        if det_levels:
            df = df[df["detection_level"].isin(det_levels)]
        src_col = f"{lvl}_region"
        tgt_col = f"{lvl}_contig"
        pairs = df[[src_col, tgt_col]].dropna().copy()
        pairs = pairs[pairs[src_col] != pairs[tgt_col]]
        pairs.columns = ["source", "target"]
        source_label = f"Region ({TAX_LEVEL_NAMES[lvl]})"
        target_label = f"Contig ({TAX_LEVEL_NAMES[lvl]})"
    else:
        df = load_hgt2(data_dir, assembly)
        src_col = f"{lvl}_contig1"
        tgt_col = f"{lvl}_contig2"
        pairs = df[[src_col, tgt_col]].dropna().copy()
        pairs = pairs[pairs[src_col] != pairs[tgt_col]]
        pairs.columns = ["source", "target"]
        source_label = f"Contig 1 ({TAX_LEVEL_NAMES[lvl]})"
        target_label = f"Contig 2 ({TAX_LEVEL_NAMES[lvl]})"

counts = (
    pairs.groupby(["source", "target"])
    .size()
    .reset_index(name="count")
)
counts = counts[counts["count"] >= min_count].sort_values("count", ascending=False)

if counts.empty:
    st.warning(
        f"No transfer pairs with count ≥ {min_count} at {TAX_LEVEL_NAMES[lvl]} level. "
        "Try lowering the minimum count or choosing a coarser taxonomic level."
    )
    st.stop()

st.caption(
    f"Showing **{len(counts):,}** transfer pairs representing **{counts['count'].sum():,}** total events "
    f"(filtered from {len(pairs):,} raw pairs)"
)

# ── Sankey ────────────────────────────────────────────────────────────────────
sources = counts["source"].tolist()
targets = counts["target"].tolist()

# Prefix source/target labels so same taxon name can appear on both sides without merging nodes
source_nodes = [f"{source_label}: {s}" for s in counts["source"].unique()]
target_nodes = [f"{target_label}: {t}" for t in counts["target"].unique()]

all_nodes = source_nodes + [n for n in target_nodes if n not in source_nodes]
node_map = {n: i for i, n in enumerate(all_nodes)}

src_idx = [node_map[f"{source_label}: {s}"] for s in sources]
tgt_idx = [node_map[f"{target_label}: {t}"] for t in targets]

# Colour nodes: sources in blue family, targets in orange family
n_src = len(source_nodes)
node_colors = (
    ["rgba(59,130,246,0.8)"] * n_src
    + ["rgba(249,115,22,0.8)"] * (len(all_nodes) - n_src)
)

fig = go.Figure(
    go.Sankey(
        arrangement="snap",
        node=dict(
            label=[n.split(": ", 1)[1] for n in all_nodes],
            pad=12,
            thickness=18,
            color=node_colors,
            hovertemplate="%{label}<extra></extra>",
        ),
        link=dict(
            source=src_idx,
            target=tgt_idx,
            value=counts["count"].tolist(),
            hovertemplate="%{source.label} → %{target.label}: %{value} events<extra></extra>",
            color="rgba(180,180,180,0.3)",
        ),
    )
)
fig.update_layout(
    title_text=f"{TAX_LEVEL_NAMES[lvl]}-level gene transfers  |  {method}",
    height=650,
    font_size=12,
)
st.plotly_chart(fig, use_container_width=True)

# ── Heatmap view (same data) ───────────────────────────────────────────────────
with st.expander("Heatmap view"):
    top_src = counts.groupby("source")["count"].sum().nlargest(30).index.tolist()
    top_tgt = counts.groupby("target")["count"].sum().nlargest(30).index.tolist()
    matrix_df = (
        counts[counts["source"].isin(top_src) & counts["target"].isin(top_tgt)]
        .pivot(index="source", columns="target", values="count")
        .fillna(0)
        .reindex(index=top_src, columns=top_tgt, fill_value=0)
    )
    fig2 = px.imshow(
        matrix_df,
        labels={"x": target_label, "y": source_label, "color": "Count"},
        color_continuous_scale="Blues",
        aspect="auto",
    )
    fig2.update_layout(height=max(400, 18 * len(top_src)))
    st.plotly_chart(fig2, use_container_width=True)

# ── Table ──────────────────────────────────────────────────────────────────────
with st.expander("Transfer pair table"):
    st.dataframe(counts.rename(columns={"source": source_label, "target": target_label}),
                 use_container_width=True, height=400)
