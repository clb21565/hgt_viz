import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import streamlit as st
import plotly.express as px
import pandas as pd
from data import load_hgt1, load_hgt2, load_dr, load_report
from config import TAX_LEVELS, TAX_LEVEL_NAMES

st.set_page_config(page_title="Overview", layout="wide")
st.title("Overview")

with st.spinner("Loading data..."):
    hgt1 = load_hgt1()
    hgt2 = load_hgt2()
    dr = load_dr()
    report = load_report()

# ── Key metrics ───────────────────────────────────────────────────────────────
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("HGT1 Regions", f"{len(hgt1):,}")
c2.metric("HGT2 ORF Pairs", f"{len(hgt2):,}")
c3.metric("Donor–Recipient Pairs", f"{len(dr):,}")
c4.metric("HGT Groups", f"{len(report):,}")
c5.metric("Contigs w/ HGT1", f"{hgt1['contig'].nunique():,}")

st.divider()

# ── Row 1: level bar + region type pie ────────────────────────────────────────
col_l, col_r = st.columns(2)

with col_l:
    st.subheader("HGT1 Events by Detection Level")
    level_counts = (
        hgt1["detection_level"]
        .value_counts()
        .reindex(TAX_LEVELS + ["unknown"], fill_value=0)
        .reset_index()
    )
    level_counts.columns = ["level", "count"]
    level_counts["Level"] = level_counts["level"].map(
        lambda x: TAX_LEVEL_NAMES.get(x, x)
    )
    fig = px.bar(
        level_counts[level_counts["count"] > 0],
        x="Level",
        y="count",
        labels={"count": "Number of Regions"},
        color="count",
        color_continuous_scale="Blues",
        text="count",
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(coloraxis_showscale=False, height=380, showlegend=False)
    st.plotly_chart(fig, use_container_width=True)

with col_r:
    st.subheader("HGT1 Region Type Distribution")
    type_counts = hgt1["region_type"].value_counts().reset_index()
    type_counts.columns = ["region_type", "count"]
    fig = px.pie(
        type_counts,
        names="region_type",
        values="count",
        color_discrete_sequence=px.colors.qualitative.Set2,
        hole=0.4,
    )
    fig.update_layout(height=380)
    st.plotly_chart(fig, use_container_width=True)

st.divider()

# ── Row 2: region length histogram + group scatter ───────────────────────────
col_l2, col_r2 = st.columns(2)

with col_l2:
    st.subheader("HGT1 Region Length Distribution")
    log_y = st.checkbox("Log y-axis", value=True, key="log_len")
    fig = px.histogram(
        hgt1,
        x="region_len",
        nbins=80,
        log_y=log_y,
        labels={"region_len": "Region Length (bp)", "count": "Count"},
        color_discrete_sequence=["#3b82f6"],
    )
    fig.update_layout(height=380, bargap=0.05)
    st.plotly_chart(fig, use_container_width=True)

with col_r2:
    st.subheader("HGT Groups: Size vs Contigs")
    color_col = st.selectbox(
        "Color by",
        ["num_distinct_phyla", "num_distinct_orders", "num_distinct_genera"],
        key="grp_color",
    )
    fig = px.scatter(
        report,
        x="num_contigs",
        y="hgt_group_size",
        color=color_col,
        hover_data=["hgt_group_id", "total_region_len", "avg_region_len", "max_region_len"],
        labels={
            "num_contigs": "Number of Contigs",
            "hgt_group_size": "Group Size (# regions)",
            color_col: color_col.replace("_", " ").title(),
        },
        color_continuous_scale="Turbo",
    )
    fig.update_layout(height=380)
    st.plotly_chart(fig, use_container_width=True)

st.divider()

# ── HGT2 contig-pair counts by phylum ────────────────────────────────────────
st.subheader("HGT2 Shared-ORF Pairs: Phylum Distribution")
st.caption("Count of shared-ORF pairs where at least one contig has the given phylum assignment.")

col_l3, col_r3 = st.columns(2)
with col_l3:
    ph1 = hgt2["p_contig1"].value_counts().head(20).reset_index()
    ph1.columns = ["Phylum", "count"]
    ph1["side"] = "Contig 1"
    ph2 = hgt2["p_contig2"].value_counts().head(20).reset_index()
    ph2.columns = ["Phylum", "count"]
    ph2["side"] = "Contig 2"
    combined = pd.concat([ph1, ph2])
    fig = px.bar(
        combined,
        x="count",
        y="Phylum",
        color="side",
        barmode="group",
        orientation="h",
        labels={"count": "Pair Count", "side": ""},
        color_discrete_sequence=["#3b82f6", "#f97316"],
    )
    fig.update_layout(height=420, yaxis={"categoryorder": "total ascending"})
    st.plotly_chart(fig, use_container_width=True)

with col_r3:
    st.subheader("")
    st.markdown(
        """
**Reading this chart:**

Each bar shows how many HGT2 (shared-ORF) pairs involve a contig taxonomically assigned to that phylum.
A pair is counted once for each contig; pairs between two contigs of the same phylum would appear in both bars (these are excluded from the Transfer Flows sankey since the taxonomy is the same).

Phyla with very high counts on *both sides* are common participants in shared-ORF HGT.
        """
    )
