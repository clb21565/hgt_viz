import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import streamlit as st
import plotly.graph_objects as go
import pandas as pd
from data import load_hgt1, load_hgt2, load_groups, dataset_sidebar
from config import TAX_LEVELS, TAX_LEVEL_NAMES
from features import Feature

st.set_page_config(page_title="Contig View", layout="wide")
st.title("Contig View")

data_dir, assembly = dataset_sidebar()

# ── Constants ──────────────────────────────────────────────────────────────────
MIN_VIS_FRAC = 0.003   # min feature width as fraction of contig length
FEAT_HALF = 0.22       # half-height of feature rectangles in y-units
BB_HALF = 0.07         # half-height of backbone

PALETTE = [
    "#e63946", "#457b9d", "#2a9d8f", "#e9c46a", "#f4a261", "#264653",
    "#8338ec", "#fb5607", "#3a86ff", "#06d6a0", "#ef476f", "#ffd166",
    "#118ab2", "#6a4c93", "#ff595e", "#ffca3a", "#8ac926", "#1982c4",
    "#c77dff", "#073b4c",
]

# ── Data loading ───────────────────────────────────────────────────────────────
@st.cache_data(show_spinner=False)
def _load_all(data_dir: str, assembly: str):
    h1 = load_hgt1(data_dir, assembly)
    h2 = load_hgt2(data_dir, assembly)
    grp = load_groups(data_dir, assembly)
    # region_idx → coordinates lookup (used for HGT2 ribbon drawing)
    # Start with groups file (covers HGT2-only contigs), then override with hgt1 (more complete)
    ridx = (
        grp.drop_duplicates("region_idx")
        .set_index("region_idx")[["contig", "start", "stop", "contig_length"]]
        .to_dict("index")
    )
    ridx.update(
        h1.drop_duplicates("region_idx")
        .set_index("region_idx")[["contig", "start", "stop", "contig_length"]]
        .to_dict("index")
    )
    return h1, h2, grp, ridx

with st.spinner("Loading data…"):
    hgt1, hgt2, groups, region_idx_lookup = _load_all(data_dir, assembly)

# ── Sidebar ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Controls")
    entry_mode = st.radio("Browse by", ["Contig", "HGT Group"])
    tax_level = st.selectbox(
        "Color regions by",
        TAX_LEVELS,
        format_func=lambda x: TAX_LEVEL_NAMES[x],
        index=1,
    )

tax_col = f"{tax_level}_region"

# ── Helper: color map ──────────────────────────────────────────────────────────
def build_color_map(series: pd.Series) -> dict:
    unique = sorted(t for t in series.dropna().unique() if str(t).strip())
    return {t: PALETTE[i % len(PALETTE)] for i, t in enumerate(unique)}

# ── Helper: min-visible expansion ─────────────────────────────────────────────
def expand_to_min_vis(start: int, end: int, contig_len: int) -> tuple[int, int, bool]:
    length = end - start
    min_len = max(1, int(contig_len * MIN_VIS_FRAC))
    if length >= min_len:
        return start, end, False
    mid = (start + end) / 2
    return int(mid - min_len / 2), int(mid + min_len / 2), True

# ── Helper: dataframe → Feature list ──────────────────────────────────────────
def df_to_features(df: pd.DataFrame, tax_col: str, color_map: dict, tax_level: str) -> list[Feature]:
    feats = []
    for _, row in df.iterrows():
        c_len = int(row["contig_length"])
        rs, re = int(row["start"]), int(row["stop"])
        vs, ve, expanded = expand_to_min_vis(rs, re, c_len)
        taxon = str(row.get(tax_col) or "").strip()
        color = color_map.get(taxon, "#999999")
        warn = " ⚠ expanded for visibility" if expanded else ""
        tt = (
            f"<b>{row.get('region_title', '')}</b><br>"
            f"Position: {rs:,}–{re:,} bp<br>"
            f"Length: {re - rs:,} bp{warn}<br>"
            f"Type: {row.get('region_type', '')}<br>"
            f"Region taxon ({TAX_LEVEL_NAMES.get(tax_level, tax_level)}): <b>{taxon or '—'}</b><br>"
            f"Contig taxon: {row.get('taxa_contig', '')}"
        )
        feats.append(Feature(
            start=vs, end=ve,
            feature_type="hgt1_region",
            label=str(row.get("region_title", "")),
            color=color, tooltip=tt,
            metadata={
                "taxon": taxon,
                "region_idx": str(row.get("region_idx", "")),
                "raw_start": rs,
                "raw_end": re,
                "expanded": expanded,
            },
        ))
    return feats

# ── Helper: draw one contig track onto a figure ───────────────────────────────
def _add_track(
    fig: go.Figure,
    features: list[Feature],
    contig_len: int,
    y_ctr: float,
    x_norm: float = 1.0,   # divide all x-coords by this (for normalized paired view)
):
    fig.add_shape(
        type="rect",
        x0=0, x1=contig_len / x_norm,
        y0=y_ctr - BB_HALF, y1=y_ctr + BB_HALF,
        fillcolor="#d1d5db", line_color="#9ca3af", line_width=0.5,
        layer="below",
    )
    for feat in features:
        x0 = feat.start / x_norm
        x1 = feat.end / x_norm
        fig.add_shape(
            type="rect",
            x0=x0, x1=x1,
            y0=y_ctr - FEAT_HALF, y1=y_ctr + FEAT_HALF,
            fillcolor=feat.color,
            line_color="rgba(0,0,0,0.15)", line_width=0.5,
        )
        # invisible scatter for hover; one point per feature at center
        fig.add_trace(go.Scatter(
            x=[(x0 + x1) / 2], y=[y_ctr],
            mode="markers",
            marker=dict(size=1, opacity=0),
            hovertemplate=feat.tooltip + "<extra></extra>",
            showlegend=False,
        ))

def _add_legend(fig: go.Figure, color_map: dict, shown_taxa: set):
    for taxon, color in color_map.items():
        if taxon in shown_taxa:
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(size=12, color=color, symbol="square"),
                name=taxon, showlegend=True,
            ))

# ── Figure builders ────────────────────────────────────────────────────────────
def make_overview(
    contig: str,
    features: list[Feature],
    contig_len: int,
    color_map: dict,
    x_range: tuple[int, int] | None = None,
    title_suffix: str = "",
) -> go.Figure:
    fig = go.Figure()
    _add_track(fig, features, contig_len, y_ctr=0.5)
    _add_legend(fig, color_map, {f.metadata["taxon"] for f in features if f.metadata.get("taxon")})
    x_range = x_range or (0, contig_len)
    fig.update_layout(
        height=160,
        margin=dict(l=0, r=10, t=32, b=32),
        xaxis=dict(range=list(x_range), title="Position (bp)", tickformat=",d"),
        yaxis=dict(range=[0, 1], visible=False),
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(orientation="v", x=1.01, y=0.5, xanchor="left", font_size=11),
        title=dict(
            text=f"{contig}  •  {contig_len:,} bp{title_suffix}",
            font_size=12, x=0,
        ),
        hovermode="closest",
    )
    return fig


def make_paired(
    contig1: str,
    contig2: str,
    feats1: list[Feature],
    feats2: list[Feature],
    len1: int,
    len2: int,
    color_map: dict,
    hgt2_pairs: pd.DataFrame,
    taxon1: str = "",
    taxon2: str = "",
) -> go.Figure:
    C1_Y, C2_Y = 2.5, 0.5

    fig = go.Figure()
    _add_track(fig, feats1, len1, y_ctr=C1_Y, x_norm=len1)
    _add_track(fig, feats2, len2, y_ctr=C2_Y, x_norm=len2)

    # Ribbons: one per unique (region_idx_contig1, region_idx_contig2) pair
    drawn: set[tuple] = set()
    for _, row in hgt2_pairs.iterrows():
        r1 = str(row["region_idx_contig1"])
        r2 = str(row["region_idx_contig2"])
        # Normalise so r1 always belongs to contig1
        if row["contig1"] == contig2:
            r1, r2 = r2, r1
        key = (r1, r2)
        if key in drawn:
            continue
        drawn.add(key)
        info1 = region_idx_lookup.get(r1)
        info2 = region_idx_lookup.get(r2)
        if not info1 or not info2:
            continue
        x1l = info1["start"] / len1
        x1r = info1["stop"] / len1
        x2l = info2["start"] / len2
        x2r = info2["stop"] / len2
        y_top = C1_Y - FEAT_HALF   # 2.28
        y_bot = C2_Y + FEAT_HALF   # 0.72
        fig.add_trace(go.Scatter(
            x=[x1l, x1r, x2r, x2l, x1l],
            y=[y_top, y_top, y_bot, y_bot, y_top],
            fill="toself",
            fillcolor="rgba(130,160,210,0.22)",
            line=dict(color="rgba(130,160,210,0.55)", width=0.8),
            hovertemplate=f"Shared ORF region<br>{r1} ↔ {r2}<extra></extra>",
            showlegend=False,
        ))

    _add_legend(
        fig, color_map,
        {f.metadata["taxon"] for f in feats1 + feats2 if f.metadata.get("taxon")},
    )
    fig.update_layout(
        height=320,
        margin=dict(l=140, r=40, t=40, b=40),
        xaxis=dict(
            range=[0, 1], title="Relative position",
            tickformat=".0%",
            tickvals=[0, 0.25, 0.5, 0.75, 1.0],
        ),
        yaxis=dict(range=[0, 3.0], visible=False),
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(orientation="v", x=1.01, y=0.5, xanchor="left", font_size=11),
        title=dict(
            text=(
                f"HGT2 Synteny: {contig1}"
                + (f" <span style='color:#6b7280;font-size:10px'>({taxon1})</span>" if taxon1 else "")
                + f" ↔ {contig2}"
                + (f" <span style='color:#6b7280;font-size:10px'>({taxon2})</span>" if taxon2 else "")
            ),
            font_size=12,
        ),
        hovermode="closest",
        annotations=[
            dict(x=-0.01, y=C1_Y + 0.15, xref="x", yref="y",
                 text=contig1, font_size=9, showarrow=False, xanchor="right"),
            dict(x=-0.01, y=C1_Y - 0.15, xref="x", yref="y",
                 text=f"<i>{taxon1 or '—'}</i>",
                 font=dict(size=9, color="#6b7280"), showarrow=False, xanchor="right"),
            dict(x=-0.01, y=C2_Y + 0.15, xref="x", yref="y",
                 text=contig2, font_size=9, showarrow=False, xanchor="right"),
            dict(x=-0.01, y=C2_Y - 0.15, xref="x", yref="y",
                 text=f"<i>{taxon2 or '—'}</i>",
                 font=dict(size=9, color="#6b7280"), showarrow=False, xanchor="right"),
            dict(x=1.01, y=C1_Y, xref="x", yref="y", text=f"{len1:,} bp",
                 font=dict(size=9, color="gray"), showarrow=False, xanchor="left"),
            dict(x=1.01, y=C2_Y, xref="x", yref="y", text=f"{len2:,} bp",
                 font=dict(size=9, color="gray"), showarrow=False, xanchor="left"),
        ],
    )
    return fig

# ── Entry selection ────────────────────────────────────────────────────────────
if entry_mode == "Contig":
    contig_list = sorted(hgt1["contig"].unique().tolist())
    selected_contig = st.selectbox(
        "Select contig",
        contig_list,
        help=f"{len(contig_list):,} contigs with HGT1 regions",
    )
    selected_contigs = [selected_contig]
else:
    group_ids = sorted(groups["hgt_group_id"].dropna().unique().tolist())
    selected_group = st.selectbox(
        "Select HGT group",
        group_ids,
        format_func=lambda x: f"Group {int(x)}",
    )
    group_rows = groups[groups["hgt_group_id"] == selected_group]
    top_contigs = (
        group_rows["contig"].value_counts().head(8).index.tolist()
    )
    selected_contigs = top_contigs
    total_contigs = group_rows["contig"].nunique()
    st.caption(
        f"Group {int(selected_group)}: {len(group_rows):,} regions across "
        f"{total_contigs:,} contigs — showing top {len(top_contigs)} by region count"
    )

# ── Per-contig rendering ───────────────────────────────────────────────────────
for contig in selected_contigs:
    regions_df = hgt1[hgt1["contig"] == contig]
    if regions_df.empty:
        st.warning(f"No HGT1 regions found for {contig}")
        continue

    contig_len = int(regions_df["contig_length"].iloc[0])
    color_map = build_color_map(
        regions_df[tax_col] if tax_col in regions_df.columns else pd.Series(dtype=str)
    )
    features = df_to_features(regions_df, tax_col, color_map, tax_level)

    # Contig background taxonomy (consistent across all rows for a given contig)
    ctax_col = f"{tax_level}_contig"
    contig_taxon = ""
    if ctax_col in regions_df.columns:
        vals = regions_df[ctax_col].dropna()
        if not vals.empty:
            contig_taxon = str(vals.mode().iloc[0])

    st.subheader(contig, divider="gray")
    st.caption(
        f"Contig background ({TAX_LEVEL_NAMES.get(tax_level, tax_level)}): "
        f"**{contig_taxon or '—'}**  •  {contig_len:,} bp"
    )

    # Overview
    fig_ov = make_overview(contig, features, contig_len, color_map)
    st.plotly_chart(fig_ov, use_container_width=True, key=f"ov_{contig}")

    # Zoom panel
    region_labels = [f.label for f in features]
    zoom_choice = st.selectbox(
        "Zoom into region",
        ["— select —"] + region_labels,
        key=f"zoom_{contig}",
        label_visibility="collapsed",
    )
    if zoom_choice != "— select —":
        feat = next((f for f in features if f.label == zoom_choice), None)
        if feat:
            rs, re = feat.metadata["raw_start"], feat.metadata["raw_end"]
            pad = max(500, int((re - rs) * 0.2))
            fig_zoom = make_overview(
                contig, features, contig_len, color_map,
                x_range=(max(0, rs - pad), min(contig_len, re + pad)),
                title_suffix=f"  |  zoom: {rs:,}–{re:,} bp",
            )
            st.plotly_chart(fig_zoom, use_container_width=True, key=f"zfig_{contig}")

    # HGT2 partners
    hgt2_hits = hgt2[(hgt2["contig1"] == contig) | (hgt2["contig2"] == contig)]
    if not hgt2_hits.empty:
        partners = (
            pd.concat([
                hgt2_hits[hgt2_hits["contig1"] == contig]["contig2"],
                hgt2_hits[hgt2_hits["contig2"] == contig]["contig1"],
            ])
            .value_counts()
            .head(20)
        )
        with st.expander(f"HGT2 partners  ({len(partners)} shown)"):
            partner_choice = st.selectbox(
                "Compare with partner",
                ["— select —"] + partners.index.tolist(),
                key=f"partner_{contig}",
            )
            if partner_choice != "— select —":
                partner_df = hgt1[hgt1["contig"] == partner_choice]

                # Fall back to groups data for HGT2-only partners (no HGT1 regions)
                if partner_df.empty:
                    grp_partner = groups[groups["contig"] == partner_choice].copy()
                    if grp_partner.empty:
                        st.warning(f"{partner_choice} has no region data available.")
                    else:
                        # Map group taxonomy columns to the _region suffix expected by df_to_features
                        for lvl in ["d", "p", "c", "o", "f", "g", "s"]:
                            grp_partner[f"{lvl}_region"] = grp_partner[lvl]
                        grp_partner["taxa_region"] = grp_partner[tax_level].fillna("")
                        grp_partner["taxa_contig"] = "(HGT group)"
                        partner_df = grp_partner
                        st.caption(f"ℹ {partner_choice} has no HGT1 regions — showing group regions instead.")

                if not partner_df.empty:
                    partner_len = int(partner_df["contig_length"].iloc[0])

                    # Partner contig background taxon
                    partner_taxon = ""
                    if ctax_col in partner_df.columns:
                        pvals = partner_df[ctax_col].dropna()
                        if not pvals.empty:
                            partner_taxon = str(pvals.mode().iloc[0])

                    combined_taxa = pd.concat([
                        regions_df[tax_col] if tax_col in regions_df.columns else pd.Series(dtype=str),
                        partner_df[tax_col] if tax_col in partner_df.columns else pd.Series(dtype=str),
                    ])
                    cmap = build_color_map(combined_taxa)
                    f1 = df_to_features(regions_df, tax_col, cmap, tax_level)
                    f2 = df_to_features(partner_df, tax_col, cmap, tax_level)
                    pairs = hgt2_hits[
                        (hgt2_hits["contig2"] == partner_choice) |
                        (hgt2_hits["contig1"] == partner_choice)
                    ]
                    fig_pair = make_paired(
                        contig, partner_choice,
                        f1, f2,
                        contig_len, partner_len,
                        cmap, pairs,
                        taxon1=contig_taxon,
                        taxon2=partner_taxon,
                    )
                    st.plotly_chart(fig_pair, use_container_width=True,
                                    key=f"pair_{contig}_{partner_choice}")
