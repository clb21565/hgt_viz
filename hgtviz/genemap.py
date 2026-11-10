"""Linear gene maps (Plotly). Genes are drawn as arrow polygons, one trace per
category, so a contig with thousands of genes stays a handful of traces."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from . import RANK_NAMES

# Categorical slots (fixed order, validated reference palette)
CARD = "#2a78d6"        # slot 1 blue
RIP = "#eb6834"         # slot 2 orange
HGT1 = "#1baf7a"        # slot 3 aqua: region band for HGT1 calls
GENE = "#a3a29c"        # neutral: unannotated gene
REGION = "#a3a29c"      # neutral: kairos region without an HGT1 call
LINK = "rgba(42,120,214,0.22)"
SEQ_BLUE = [[0, "#cde2fb"], [0.25, "#86b6ef"], [0.5, "#3987e5"], [0.75, "#1c5cab"], [1, "#0d366b"]]

GENE_H = 0.28      # half-height of a gene arrow
BAND_H = 0.42      # half-height of a region band


def gene_category(genes: pd.DataFrame, ann: pd.DataFrame | None) -> pd.DataFrame:
    """Attach annotation columns and a drawing category to each gene."""
    g = genes.copy()
    g["category"] = "gene"
    g["annotation"] = ""
    if ann is not None and len(ann):
        a = ann[ann["orf_id"].isin(g["orf_id"])]
        for src, cat in [("ripdb", "RIP-DB"), ("card", "CARD")]:  # CARD wins if both
            s = a[a["source"] == src].drop_duplicates("orf_id").set_index("orf_id")
            hit = g["orf_id"].isin(s.index)
            g.loc[hit, "category"] = cat
            txt = s["label"] + " (" + s["pident"].round(1).astype(str) + "% id"
            if src == "card":
                txt = txt + np.where(s["category"].notna(), ", " + s["category"].fillna(""), "")
            g.loc[hit, "annotation"] = g.loc[hit, "orf_id"].map(txt + ")")
    return g


def _arrows(df, y, xf, head):
    """Polygon coordinates for all genes in df, separated by None."""
    xs, ys = [], []
    for s, e, st in zip(df["start"], df["end"], df["strand"]):
        a, b = sorted((xf(s), xf(e)))
        st = st * (1 if xf(e) >= xf(s) else -1)
        h = min(head, (b - a) * 0.5)
        if st >= 0:
            px = [a, b - h, b, b - h, a, a]
        else:
            px = [b, a + h, a, a + h, b, b]
        py = [y - GENE_H, y - GENE_H, y, y + GENE_H, y + GENE_H, y - GENE_H]
        xs += px + [None]
        ys += py + [None]
    return xs, ys


def add_track(fig, genes, regions, y, label, length, xf=lambda x: x, x_range=None):
    """Draw one contig: backbone, region bands, genes, gene hover points."""
    lo, hi = sorted((xf(0), xf(length)))
    span = (x_range[1] - x_range[0]) if x_range else (hi - lo)
    head = span * 0.006
    fig.add_shape(type="line", x0=lo, x1=hi, y0=y, y1=y,
                  line=dict(color="#a3a29c", width=1.5), layer="below")
    for _, r in regions.iterrows():
        a, b = sorted((xf(r["start"]), xf(r["end"])))
        hgt1 = pd.notna(r.get("hgt1_rank"))
        fig.add_shape(type="rect", x0=a, x1=b, y0=y - BAND_H, y1=y + BAND_H,
                      fillcolor=HGT1 if hgt1 else REGION, opacity=0.22, line_width=0,
                      layer="below")
    # region hover: invisible points along the band
    if len(regions):
        mids = [(xf(s) + xf(e)) / 2 for s, e in zip(regions["start"], regions["end"])]
        txt = [_region_hover(r) for _, r in regions.iterrows()]
        fig.add_trace(go.Scatter(x=mids, y=[y + BAND_H] * len(mids), mode="markers",
                                 marker=dict(size=14, opacity=0, color=HGT1),
                                 hovertemplate="%{text}<extra></extra>", text=txt,
                                 showlegend=False))
    if genes is not None and len(genes):
        for cat, color in [("gene", GENE), ("RIP-DB", RIP), ("CARD", CARD)]:
            g = genes[genes["category"] == cat]
            if g.empty:
                continue
            xs, ys = _arrows(g, y, xf, head)
            fig.add_trace(go.Scatter(x=xs, y=ys, fill="toself", mode="lines",
                                     fillcolor=color, line=dict(color=color, width=0.5),
                                     name=cat, legendgroup=cat, hoverinfo="skip",
                                     showlegend=False))
        mid = [(xf(s) + xf(e)) / 2 for s, e in zip(genes["start"], genes["end"])]
        text = [f"<b>{o}</b><br>{s + 1:,}–{e:,} ({'+' if st > 0 else '−'})"
                + (f"<br>{an}" if an else "")
                + (f"<br>region {ri}" if isinstance(ri, str) else "")
                for o, s, e, st, an, ri in zip(genes["orf_id"], genes["start"], genes["end"],
                                               genes["strand"], genes["annotation"],
                                               genes.get("region_idx", [None] * len(genes)))]
        fig.add_trace(go.Scatter(x=mid, y=[y] * len(mid), mode="markers",
                                 marker=dict(size=10, opacity=0),
                                 hovertemplate="%{text}<extra></extra>", text=text,
                                 showlegend=False))
    fig.add_annotation(x=0, xref="paper", y=y + BAND_H + 0.04, text=label, showarrow=False,
                       xanchor="left", yanchor="bottom", font=dict(size=12))


def _region_hover(r):
    rank = r.get("hgt1_rank")
    lines = [f"<b>{r['region_idx']}</b>", f"{r['start'] + 1:,}–{r['end']:,} bp ({r['region_type']})"]
    if pd.notna(rank):
        lines.append(f"HGT1 at {RANK_NAMES.get(rank, rank).lower()} level")
    if pd.notna(r.get("region_taxon")):
        lines.append(f"region taxon: {r['region_taxon']}")
    if pd.notna(r.get("group_id")):
        lines.append(f"group {int(r['group_id'])} ({int(r['group_size'])} regions)")
    return "<br>".join(lines)


def _legend(fig, cats):
    for name, color in [("CARD hit", CARD), ("RIP-DB hit", RIP), ("other gene", GENE),
                        ("HGT1 region", HGT1), ("other kairos region", REGION)]:
        if name in cats:
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=name,
                                     marker=dict(size=11, color=color, symbol="square",
                                                 opacity=0.35 if "region" in name else 1)))


def _layout(fig, height, x_range, xtitle="Position (bp)"):
    fig.update_layout(
        height=height, margin=dict(l=10, r=10, t=10, b=40),
        xaxis=dict(title=xtitle, tickformat=",d", range=x_range, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, fixedrange=True),
        legend=dict(orientation="h", y=-0.25, x=0),
        hovermode="closest", dragmode="zoom",
        plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
    )


def contig_map(contig, length, genes, regions, links=None, x_range=None):
    """One contig, with shared segments to partner contigs as ticks underneath."""
    fig = go.Figure()
    add_track(fig, genes, regions, 1.0, f"{contig} · {length:,} bp", length, x_range=x_range)
    if links is not None and len(links):
        # one tick per distinct shared segment, stacked to avoid overlap
        segs = (links.groupby(["start", "end"])
                .agg(n=("partner", "nunique"), partners=("partner", lambda s: ", ".join(sorted(set(s))[:6])))
                .reset_index())
        xs, ys = [], []
        for i, (_, s) in enumerate(segs.iterrows()):
            yy = 0.35 - 0.12 * (i % 3)
            xs += [s["start"], s["end"], None]
            ys += [yy, yy, None]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=CARD, width=4),
                                 hoverinfo="skip", showlegend=False))
        fig.add_trace(go.Scatter(
            x=(segs["start"] + segs["end"]) / 2,
            y=[0.35 - 0.12 * (i % 3) for i in range(len(segs))], mode="markers",
            marker=dict(size=12, opacity=0),
            text=[f"shared with {n} contig(s)<br>{p}" for n, p in zip(segs["n"], segs["partners"])],
            hovertemplate="%{text}<extra></extra>", showlegend=False))
        fig.add_annotation(x=0, y=0.5, text="shared segments", showarrow=False, xanchor="left",
                           font=dict(size=11, color="#52514e"))
    cats = {"other gene"}
    if genes is not None:
        cats |= {"CARD hit" if c == "CARD" else "RIP-DB hit" for c in genes["category"].unique()
                 if c != "gene"}
    if len(regions):
        cats |= {"HGT1 region" if pd.notna(x) else "other kairos region" for x in regions["hgt1_rank"]}
    _legend(fig, cats)
    _layout(fig, 260, x_range or [0, length])
    fig.update_yaxes(range=[0, 1.65])
    return fig


def paired_map(c1, len1, genes1, regions1, c2, len2, genes2, regions2, link, pad=5000):
    """Two contigs aligned on their shared segment. The partner is shifted (and
    flipped when the copies are in opposite orientation) so the segment lines up."""
    s1, e1 = int(link["start"]), int(link["end"])
    s2, e2 = int(link["partner_start"]), int(link["partner_end"])
    if int(link.get("strand", 1) or 1) >= 0:
        def xf2(x):
            return s1 + (x - s2)
    else:
        def xf2(x):
            return s1 + (e2 - x)
    lo = max(min(0, xf2(0), xf2(len2)), s1 - pad)
    hi = min(max(len1, xf2(0), xf2(len2)), e1 + pad)
    x_range = [lo, hi]
    fig = go.Figure()
    a, b = sorted((xf2(s2), xf2(e2)))
    y2 = -0.6
    fig.add_trace(go.Scatter(x=[s1, e1, b, a, s1], y=[1.0 - BAND_H, 1.0 - BAND_H, y2 + BAND_H,
                                                      y2 + BAND_H, 1.0 - BAND_H],
                             fill="toself", fillcolor=LINK, line=dict(width=0), mode="lines",
                             hoverinfo="skip", showlegend=False))
    add_track(fig, genes1, regions1, 1.0, f"{c1} · {len1:,} bp", len1, x_range=x_range)
    strand_note = " (reverse complement)" if int(link.get("strand", 1) or 1) < 0 else ""
    add_track(fig, genes2, regions2, y2, f"{c2} · {len2:,} bp{strand_note}", len2, xf=xf2,
              x_range=x_range)
    cats = {"other gene", "HGT1 region", "other kairos region"}
    for g in (genes1, genes2):
        if g is not None:
            cats |= {"CARD hit" if c == "CARD" else "RIP-DB hit" for c in g["category"].unique()
                     if c != "gene"}
    _legend(fig, cats)
    _layout(fig, 340, x_range, xtitle=f"Position on {c1} (bp)")
    fig.update_yaxes(range=[y2 - 0.55, 1.65])
    return fig
