import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import streamlit as st
import plotly.express as px
import pandas as pd
from data import load_report, load_groups, dataset_sidebar
from config import TAX_LEVELS, TAX_LEVEL_NAMES

st.set_page_config(page_title="HGT Groups", layout="wide")
st.title("HGT Groups")

data_dir, assembly = dataset_sidebar()
st.markdown(
    "Groups cluster HGT regions from multiple contigs that likely share a common transferred element. "
    "Explore group size, taxonomic diversity, and region characteristics."
)

with st.spinner("Loading group data..."):
    report = load_report(data_dir, assembly)

# ── Controls ──────────────────────────────────────────────────────────────────
col1, col2, col3 = st.columns(3)
with col1:
    size_metric = st.selectbox(
        "Size metric",
        ["total_region_len", "hgt_group_size", "num_regions", "num_contigs"],
        format_func=lambda x: x.replace("_", " ").title(),
    )
with col2:
    color_metric = st.selectbox(
        "Color metric",
        ["num_distinct_phyla", "num_distinct_classes", "num_distinct_orders",
         "num_distinct_genera", "num_distinct_species"],
        format_func=lambda x: x.replace("_", " ").title(),
    )
with col3:
    top_n = st.slider("Top N groups to display", 20, min(500, len(report)), 100)

top = report.nlargest(top_n, size_metric).copy()
top["group_label"] = top["hgt_group_id"].astype(str)

# ── Scatter ───────────────────────────────────────────────────────────────────
st.subheader("Group Scatter")
col_s, col_s2 = st.columns([3, 1])
with col_s2:
    x_metric = st.selectbox(
        "X axis",
        ["num_contigs", "num_regions", "avg_region_len", "max_region_len"],
        format_func=lambda x: x.replace("_", " ").title(),
    )
    log_x = st.checkbox("Log X", value=False)
    log_y = st.checkbox("Log Y", value=False)

with col_s:
    fig = px.scatter(
        top,
        x=x_metric,
        y=size_metric,
        size=size_metric,
        color=color_metric,
        hover_data={
            "group_label": True,
            "hgt_group_size": True,
            "num_regions": True,
            "num_contigs": True,
            "total_region_len": True,
            "avg_region_len": ":.0f",
            "max_region_len": True,
            "num_distinct_phyla": True,
        },
        labels={
            x_metric: x_metric.replace("_", " ").title(),
            size_metric: size_metric.replace("_", " ").title(),
            color_metric: color_metric.replace("_", " ").title(),
            "group_label": "Group ID",
        },
        color_continuous_scale="Turbo",
        log_x=log_x,
        log_y=log_y,
    )
    fig.update_layout(height=480)
    st.plotly_chart(fig, use_container_width=True)

st.divider()

# ── Treemap ───────────────────────────────────────────────────────────────────
st.subheader("Treemap")
fig = px.treemap(
    top,
    path=["group_label"],
    values=size_metric,
    color=color_metric,
    hover_data=["num_contigs", "num_regions", "avg_region_len", "max_region_len"],
    color_continuous_scale="RdYlBu_r",
    title=f"Top {top_n} HGT groups — sized by {size_metric.replace('_', ' ')}",
)
fig.update_traces(textinfo="label+value")
fig.update_layout(height=560)
st.plotly_chart(fig, use_container_width=True)

st.divider()

# ── Diversity breakdown bar ────────────────────────────────────────────────────
st.subheader("Taxonomic Diversity in Top Groups")
diversity_cols = [
    "num_distinct_phyla", "num_distinct_classes", "num_distinct_orders",
    "num_distinct_families", "num_distinct_genera", "num_distinct_species",
]
top25 = report.nlargest(25, size_metric)[["hgt_group_id"] + diversity_cols].copy()
top25["hgt_group_id"] = top25["hgt_group_id"].astype(str)
melted = top25.melt(id_vars="hgt_group_id", value_vars=diversity_cols,
                    var_name="level", value_name="distinct_taxa")
melted["level"] = melted["level"].str.replace("num_distinct_", "").str.title()

fig = px.bar(
    melted,
    x="hgt_group_id",
    y="distinct_taxa",
    color="level",
    barmode="group",
    labels={"hgt_group_id": "Group ID", "distinct_taxa": "Distinct Taxa", "level": "Level"},
    color_discrete_sequence=px.colors.qualitative.Plotly,
)
fig.update_layout(height=400, xaxis_tickangle=-45)
st.plotly_chart(fig, use_container_width=True)

# ── Table ──────────────────────────────────────────────────────────────────────
with st.expander("Full group table (sortable)"):
    display_cols = [
        "hgt_group_id", "hgt_group_size", "num_regions", "num_contigs",
        "total_region_len", "avg_region_len", "max_region_len",
        "num_distinct_phyla", "num_distinct_orders", "num_distinct_genera", "num_distinct_species",
    ]
    st.dataframe(
        report[[c for c in display_cols if c in report.columns]]
        .sort_values("hgt_group_size", ascending=False)
        .reset_index(drop=True),
        use_container_width=True,
        height=420,
    )
