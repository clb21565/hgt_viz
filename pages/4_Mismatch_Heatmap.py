import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import streamlit as st
import plotly.express as px
import pandas as pd
from data import load_hgt1, dataset_sidebar
from config import TAX_LEVELS, TAX_LEVEL_NAMES

st.set_page_config(page_title="Mismatch Heatmap", layout="wide")
st.title("Taxonomic Mismatch Heatmap")

data_dir, assembly = dataset_sidebar()
st.markdown(
    "Each cell counts how many **HGT1 regions** have a given *region taxonomy* (Y axis, the putative donor lineage) "
    "embedded in a contig with a different *contig taxonomy* (X axis, the recipient organism). "
    "Rows and columns are ordered by total event count."
)

with st.spinner("Loading HGT1 data..."):
    hgt1 = load_hgt1(data_dir, assembly)

# ── Controls ──────────────────────────────────────────────────────────────────
col1, col2, col3 = st.columns(3)
with col1:
    lvl = st.selectbox(
        "Taxonomic level to compare",
        TAX_LEVELS,
        format_func=lambda x: TAX_LEVEL_NAMES[x],
        index=1,  # Phylum default
    )
with col2:
    top_n = st.slider("Top N taxa per axis", 5, 40, 20)
with col3:
    det_filter = st.multiselect(
        "Filter by detection level",
        options=TAX_LEVELS,
        default=TAX_LEVELS,
        format_func=lambda x: TAX_LEVEL_NAMES[x],
    )

df = hgt1[hgt1["detection_level"].isin(det_filter)] if det_filter else hgt1

r_col = f"{lvl}_region"
c_col = f"{lvl}_contig"

pairs = df[[r_col, c_col]].dropna().copy()
pairs = pairs[pairs[r_col] != pairs[c_col]]
pairs.columns = ["region_taxon", "contig_taxon"]

if pairs.empty:
    st.warning("No mismatches found at this level with the current filters.")
    st.stop()

# ── Build matrix ──────────────────────────────────────────────────────────────
top_region = pairs["region_taxon"].value_counts().head(top_n).index.tolist()
top_contig = pairs["contig_taxon"].value_counts().head(top_n).index.tolist()

filtered = pairs[pairs["region_taxon"].isin(top_region) & pairs["contig_taxon"].isin(top_contig)]

matrix = (
    filtered.groupby(["region_taxon", "contig_taxon"])
    .size()
    .unstack(fill_value=0)
    .reindex(index=top_region, columns=top_contig, fill_value=0)
)

st.caption(
    f"Showing top {top_n} taxa on each axis — "
    f"**{len(filtered):,}** of {len(pairs):,} total mismatch events."
)

# ── Heatmap ───────────────────────────────────────────────────────────────────
log_scale = st.checkbox("Log color scale", value=True)

import numpy as np
z = matrix.values
if log_scale:
    z_plot = np.log1p(z)
    colorbar_title = "log(count + 1)"
else:
    z_plot = z
    colorbar_title = "Count"

fig = px.imshow(
    z_plot,
    x=matrix.columns.tolist(),
    y=matrix.index.tolist(),
    labels={
        "x": f"Contig taxon ({TAX_LEVEL_NAMES[lvl]})",
        "y": f"Region taxon ({TAX_LEVEL_NAMES[lvl]})",
        "color": colorbar_title,
    },
    color_continuous_scale="Blues",
    aspect="auto",
    text_auto=False,
)
fig.update_xaxes(tickangle=-45)
height = max(400, 22 * top_n)
fig.update_layout(height=height)
st.plotly_chart(fig, use_container_width=True)

# ── Top pairs table ───────────────────────────────────────────────────────────
st.subheader("Top Transfer Pairs")
top_pairs = (
    pairs.groupby(["region_taxon", "contig_taxon"])
    .size()
    .reset_index(name="count")
    .sort_values("count", ascending=False)
    .head(50)
    .reset_index(drop=True)
)
top_pairs.columns = [
    f"Region Taxon ({TAX_LEVEL_NAMES[lvl]})",
    f"Contig Taxon ({TAX_LEVEL_NAMES[lvl]})",
    "Count",
]
st.dataframe(top_pairs, use_container_width=True, height=350)
