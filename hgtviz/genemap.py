"""Linear gene maps (Plotly). Genes are drawn as arrow polygons, one trace per
category, so a contig with thousands of genes stays a handful of traces."""
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from . import RANK_NAMES, RANKS
from . import cog as COG

# Categorical slots (fixed order, validated reference palette)
CARD = "#2a78d6"        # slot 1 blue
RIP = "#eb6834"         # slot 2 orange
HGT1 = "#1baf7a"        # slot 3 aqua: region band for HGT1 calls
GENE = "#a3a29c"        # neutral: unannotated gene
REGION = "#a3a29c"      # neutral: kairos region without an HGT1 call
LINK = "rgba(42,120,214,0.22)"
# COG mode: categorical slots 1, 2, 4, 5 in fixed order (slot 3 aqua stays the HGT1
# region band); the two uninformative groups are neutral greys.
COG_COLORS = {COG.MOBILE: "#2a78d6", COG.INFO: "#eb6834", COG.CELL: "#eda100",
              COG.METAB: "#e87ba4", COG.POOR: "#a3a29c", COG.NONE: "#d6d5d0"}
SEQ_BLUE = [[0, "#cde2fb"], [0.25, "#86b6ef"], [0.5, "#3987e5"], [0.75, "#1c5cab"], [1, "#0d366b"]]

# ORF-cluster mode: the four categorical hues that pass every pairwise check (clusters
# sit side by side in any order), each also labelled with its number on the map.
CLUSTER_COLORS = ["#2a78d6", "#eda100", "#008300", "#e34948"]
SHARED_ORF = "#6b6a65"    # identical copy on another contig, beyond the 4 coloured clusters
UNIQUE_ORF = "#d6d5d0"    # no identical copy anywhere in the catalogue
SEGMENT = "#6b6a65"   # neutral: shared-segment ticks (blue is a gene colour in both modes)

# HGT-evidence mode: one blue ramp, darker = deeper taxonomic disagreement. Seven
# ranks do not fit as distinguishable steps (validated, --ordinal), so family, genus
# and species share the lightest step; hover and region labels give the exact rank.
HGT_LEVELS = {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family or below"}
HGT_RAMP = {"d": "#0d366b", "p": "#1c5cab", "c": "#2a78d6", "o": "#5598e7", "f": "#86b6ef"}
NO_HGT_REGION = "#a3a29c"   # in a kairos region without an HGT call
NO_HGT = "#d6d5d0"          # outside kairos regions


def _level(rank):
    return None if not isinstance(rank, str) else ("f" if rank in "fgs" else rank)


def _coarsest(ranks):
    idx = [RANKS.index(r) for r in ranks if isinstance(r, str) and r in RANKS]
    return RANKS[min(idx)] if idx else None


def _rname(r):
    return RANK_NAMES.get(r, r).lower() if isinstance(r, str) else ""


GENE_H = 0.28      # half-height of a gene arrow
BAND_H = 0.42      # half-height of a region band


def gene_category(genes: pd.DataFrame, ann: pd.DataFrame | None,
                  eggnog: pd.DataFrame | None = None) -> pd.DataFrame:
    """Attach annotation columns and drawing categories to each gene: `category`
    (CARD / RIP-DB / gene) and `cog_group`, plus hover text."""
    g = genes.copy()
    g["category"] = "gene"
    g["annotation"] = ""
    g["cog_group"] = COG.NONE
    g["eggnog"] = ""
    g["clcat"] = "unique"
    if "cl_contigs" in g.columns:
        g["cl_text"] = np.where(
            g["cluster"] >= 0,
            "ORF cluster #" + g["cluster"].astype(str) + ": " + g["cl_members"].astype(str)
            + " copies on " + g["cl_contigs"].astype(str) + " contig(s)", "")
    else:
        g["cl_text"] = ""
    if eggnog is not None and len(eggnog):
        e = eggnog.drop_duplicates("orf_id").set_index("orf_id")
        hit = g["orf_id"].isin(e.index)
        o = g.loc[hit, "orf_id"]
        g.loc[hit, "cog_group"] = o.map(e["cog_group"]).fillna(COG.NONE).to_numpy()

        def txt(r):
            name = r["preferred_name"] if isinstance(r["preferred_name"], str) else ""
            desc = r["description"] if isinstance(r["description"], str) else ""
            bits = [f"eggNOG: {name + ' · ' if name else ''}{desc}".rstrip(" ·")]
            if isinstance(r["cog_category"], str):
                bits.append("COG " + ", ".join(COG.label(c) for c in COG.letters(r["cog_category"])))
            extra = [x for x in (r["kegg_ko"], r["pfams"]) if isinstance(x, str)]
            if extra:
                bits.append(" · ".join(extra)[:120])
            return "<br>".join(bits)
        g.loc[hit, "eggnog"] = o.map(e.apply(txt, axis=1)).to_numpy()
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


def annotate_hgt(genes, regions, links, host=None):
    """HGT evidence per region and per gene, regardless of type.

    HGT1: the gene's region has an HGT1 call (region taxonomy differs from the host
    contig's). HGT2: the gene overlaps a shared segment with a contig of another
    lineage. A gene's `hgt_rank` is the coarsest rank over both; `hgt_cat` is its
    ramp level, or 'region' (in a kairos region, no call) / 'none'. Adds hover text
    (`hgt_text`) and, on regions, a short label (`hgt_label`) such as
    'HGT1 phylum · HGT2 phylum (13)'."""
    regions = regions.copy()
    L = links[links["is_hgt2"]] if links is not None and len(links) else None

    def host_at(r):
        return host[r] if host is not None and r in host and pd.notna(host[r]) else "?"

    # regions
    labels, h2r, details = [], [], []
    for _, r in regions.iterrows():
        parts, det = [], []
        r1 = r.get("hgt1_rank")
        if isinstance(r1, str):
            parts.append(f"HGT1 {_rname(r1)}")
            det.append(f"HGT1 at {_rname(r1)}: region {r.get(f'r_{r1}', '?')} vs host {host_at(r1)}")
        rr = None
        if L is not None:
            lr = L[L["region"] == r["region_idx"]]
            if len(lr):
                rr = _coarsest(lr["discordant_rank"])
                n = lr["partner"].nunique()
                parts.append(f"HGT2 {_rname(rr)} ({n})")
                det.append(f"HGT2 with {n} contig(s), coarsest at {_rname(rr)}: "
                           + _taxa(lr[lr["discordant_rank"] == rr], rr)
                           + (f"; best identity {lr['identity'].max():.1%}"
                              if "identity" in lr and lr["identity"].notna().any() else ""))
        h2r.append(rr)
        labels.append(" · ".join(parts))
        details.append("<br>".join(det))
    regions["h2_rank"] = h2r
    regions["hgt_label"] = labels
    regions["hgt_detail"] = details
    if genes is None:
        return None, regions

    # genes
    g = genes.copy()
    reg = regions.set_index("region_idx")
    h1 = g["region_idx"].map(reg["hgt1_rank"]) if "region_idx" in g else pd.Series(None, index=g.index)
    h2 = [None] * len(g)
    h2txt = [""] * len(g)
    if L is not None and {"start", "end"} <= set(L.columns) and L["start"].notna().any():
        Ls = L.dropna(subset=["start", "end"])
        S, E = Ls["start"].to_numpy(), Ls["end"].to_numpy()
        rk = Ls["discordant_rank"].map(lambda x: RANKS.index(x) if x in RANKS else 99).to_numpy()
        ov = (S[None, :] < g["end"].to_numpy()[:, None]) & (E[None, :] > g["start"].to_numpy()[:, None])
        for i in np.flatnonzero(ov.any(axis=1)):
            hit = Ls[ov[i]]
            rr = RANKS[rk[ov[i]].min()] if rk[ov[i]].min() < 99 else None
            h2[i] = rr
            h2txt[i] = (f"HGT2: shared with {hit['partner'].nunique()} contig(s) of another lineage; "
                        f"coarsest at {_rname(rr)}: {_taxa(hit[hit['discordant_rank'] == rr], rr)}")
    g["hgt_rank"] = [_coarsest([a, b]) for a, b in zip(h1, h2)]
    g["hgt_cat"] = [_level(r) or ("region" if isinstance(ri, str) else "none")
                    for r, ri in zip(g["hgt_rank"], g.get("region_idx", [None] * len(g)))]
    h1txt = [f"HGT1: region {reg.at[ri, f'r_{r}'] if f'r_{r}' in reg else '?'} vs host {host_at(r)} "
             f"(differ at {_rname(r)})" if isinstance(r, str) else ""
             for r, ri in zip(h1, g.get("region_idx", [None] * len(g)))]
    g["hgt_text"] = ["<br>".join(x for x in (a, b) if x) for a, b in zip(h1txt, h2txt)]
    return g, regions


def _taxa(df, rank, k=3):
    if df.empty or not isinstance(rank, str) or f"p_{rank}" not in df:
        return "?"
    v = df.drop_duplicates("partner")[f"p_{rank}"].value_counts()
    return ", ".join(f"{t} ({n})" for t, n in v.head(k).items()) + (" …" if len(v) > k else "")


def assign_clusters(frames, paired=False, n=len(CLUSTER_COLORS)):
    """Pick up to n ORF clusters to colour in this view and set each gene's `clcat`
    (cl1..cln / shared / unique). Single contig: clusters with copies on the most
    other contigs. Paired: clusters present on both contigs, most genes first.
    Returns {cluster id: (number, colour, cl_contigs)}."""
    frames = [f for f in frames if f is not None]
    if not frames or "cluster" not in frames[0].columns:
        return {}
    for f in frames:
        f["clcat"] = np.where(f["cl_contigs"] > 1, "shared", "unique")
    if paired and len(frames) == 2:
        both = set(frames[0]["cluster"]) & set(frames[1]["cluster"]) - {-1}
        cand = pd.concat(frames)
        cand = cand[cand["cluster"].isin(both)]
        order = (cand.groupby("cluster").agg(n=("orf_id", "size"), pos=("start", "min"))
                 .sort_values(["n", "pos"], ascending=[False, True]).index)
    else:
        f = frames[0]
        cand = f[(f["cl_contigs"] > 1)]
        order = (cand.groupby("cluster").agg(k=("cl_contigs", "max"), pos=("start", "min"))
                 .sort_values(["k", "pos"], ascending=[False, True]).index)
    pick = {}
    for i, cl in enumerate(order[:n]):
        k = int(pd.concat(frames).loc[lambda d: d["cluster"] == cl, "cl_contigs"].iloc[0])
        pick[cl] = (i + 1, CLUSTER_COLORS[i], k)
    for f in frames:
        for cl, (i, _, _) in pick.items():
            f.loc[f["cluster"] == cl, "clcat"] = f"cl{i}"
    return pick


def add_track(fig, genes, regions, y, label, length, xf=lambda x: x, x_range=None,
              color_by="hits"):
    """Draw one contig: backbone, region bands, genes, gene hover points."""
    lo, hi = sorted((xf(0), xf(length)))
    span = (x_range[1] - x_range[0]) if x_range else (hi - lo)
    head = span * 0.006
    fig.add_shape(type="line", x0=lo, x1=hi, y0=y, y1=y,
                  line=dict(color="#a3a29c", width=1.5), layer="below")
    for _, r in regions.iterrows():
        a, b = sorted((xf(r["start"]), xf(r["end"])))
        hgt1 = pd.notna(r.get("hgt1_rank"))
        band = REGION if color_by == "hgt" else (HGT1 if hgt1 else REGION)
        fig.add_shape(type="rect", x0=a, x1=b, y0=y - BAND_H, y1=y + BAND_H,
                      fillcolor=band, opacity=0.22, line_width=0, layer="below")
    if color_by == "hgt" and "hgt_label" in regions and (regions["hgt_label"] != "").any():
        rl = regions[regions["hgt_label"] != ""]
        fig.add_trace(go.Scatter(
            x=[(xf(s) + xf(e)) / 2 for s, e in zip(rl["start"], rl["end"])],
            y=[y + BAND_H + 0.1] * len(rl), mode="text", text=rl["hgt_label"],
            textfont=dict(size=11, color="#52514e"), hoverinfo="skip", showlegend=False))
    # region hover: invisible points along the band
    if len(regions):
        mids = [(xf(s) + xf(e)) / 2 for s, e in zip(regions["start"], regions["end"])]
        txt = [_region_hover(r) for _, r in regions.iterrows()]
        fig.add_trace(go.Scatter(x=mids, y=[y + BAND_H] * len(mids), mode="markers",
                                 marker=dict(size=14, opacity=0, color=HGT1),
                                 hovertemplate="%{text}<extra></extra>", text=txt,
                                 showlegend=False))
    if genes is not None and len(genes):
        if color_by == "cog":
            layers = [(genes["cog_group"] == grp, c) for grp, c in COG_COLORS.items()]
        elif color_by == "hgt":
            layers = [(genes["hgt_cat"] == "none", NO_HGT), (genes["hgt_cat"] == "region", NO_HGT_REGION)]
            layers += [(genes["hgt_cat"] == k, c) for k, c in reversed(list(HGT_RAMP.items()))]
        elif color_by == "cluster":
            layers = [(genes["clcat"] == "unique", UNIQUE_ORF), (genes["clcat"] == "shared", SHARED_ORF)]
            layers += [(genes["clcat"] == f"cl{i + 1}", c) for i, c in enumerate(CLUSTER_COLORS)]
        else:
            layers = [(genes["category"] == cat, c) for cat, c in
                      [("gene", GENE), ("RIP-DB", RIP), ("CARD", CARD)]]
        for sel, color in layers:
            g = genes[sel]
            if g.empty:
                continue
            xs, ys = _arrows(g, y, xf, head)
            fig.add_trace(go.Scatter(x=xs, y=ys, fill="toself", mode="lines",
                                     fillcolor=color, line=dict(color=color, width=0.5),
                                     hoverinfo="skip",
                                     showlegend=False))
        if color_by == "cluster":
            lab = genes[genes["clcat"].str.startswith("cl")]
            if len(lab):
                fig.add_trace(go.Scatter(
                    x=[(xf(s) + xf(e)) / 2 for s, e in zip(lab["start"], lab["end"])],
                    y=[y + GENE_H + 0.1] * len(lab), mode="text", text=lab["clcat"].str[2:],
                    textfont=dict(size=11, color="#52514e"), hoverinfo="skip", showlegend=False))
        mid = [(xf(s) + xf(e)) / 2 for s, e in zip(genes["start"], genes["end"])]
        text = [f"<b>{o}</b><br>{s + 1:,}–{e:,} ({'+' if st > 0 else '−'})"
                + (f"<br>{hg}" if hg else "")
                + (f"<br>{an}" if an else "")
                + (f"<br>{eg}" if eg else "")
                + (f"<br>{cl}" if cl else "")
                + (f"<br>region {ri}" if isinstance(ri, str) else "")
                for o, s, e, st, hg, an, eg, cl, ri in zip(
                    genes["orf_id"], genes["start"], genes["end"], genes["strand"],
                    genes.get("hgt_text", [""] * len(genes)),
                    genes["annotation"], genes.get("eggnog", [""] * len(genes)),
                    genes.get("cl_text", [""] * len(genes)),
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
    if isinstance(r.get("hgt_detail"), str) and r["hgt_detail"]:
        lines.append(r["hgt_detail"])
    elif pd.notna(rank):
        lines.append(f"HGT1 at {RANK_NAMES.get(rank, rank).lower()} level")
    if pd.notna(r.get("region_taxon")):
        lines.append(f"region taxon: {r['region_taxon']}")
    if pd.notna(r.get("group_id")):
        lines.append(f"group {int(r['group_id'])} ({int(r['group_size'])} regions)")
    return "<br>".join(lines)


def _legend(fig, cats, clusters=None):
    for cl, (i, color, k) in (clusters or {}).items():
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                                 name=f"{i}: ORF cluster #{cl} ({k} contigs)",
                                 marker=dict(size=11, color=color, symbol="square")))
    for name, color in ([("CARD hit", CARD), ("RIP-DB hit", RIP), ("other gene", GENE)]
                        + list(COG_COLORS.items())
                        + [(f"HGT: {n}", HGT_RAMP[k]) for k, n in HGT_LEVELS.items()]
                        + [("kairos region, no HGT call", NO_HGT_REGION), ("no HGT evidence", NO_HGT)]
                        + [("identical copy on another contig", SHARED_ORF),
                           ("unique ORF (no identical copy)", UNIQUE_ORF)]
                        + [("HGT1 region", HGT1), ("other kairos region", REGION),
                           ("kairos region", REGION)]):
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


def _gene_cats(genes, color_by):
    if genes is None:
        return set()
    if color_by == "cog":
        return set(genes["cog_group"].unique())
    if color_by == "hgt":
        return {f"HGT: {HGT_LEVELS[c]}" if c in HGT_LEVELS else
                {"region": "kairos region, no HGT call", "none": "no HGT evidence"}[c]
                for c in genes["hgt_cat"].unique()}
    if color_by == "cluster":
        names = {"unique": "unique ORF (no identical copy)",
                 "shared": "identical copy on another contig"}
        return {names[c] for c in genes["clcat"].unique() if c in names}
    return {"other gene"} | {"CARD hit" if c == "CARD" else "RIP-DB hit"
                             for c in genes["category"].unique() if c != "gene"}


def contig_map(contig, length, genes, regions, links=None, x_range=None, color_by="hits"):
    """One contig, with shared segments to partner contigs as ticks underneath."""
    fig = go.Figure()
    if genes is not None:
        genes = genes.copy()
    clusters = assign_clusters([genes]) if color_by == "cluster" else {}
    add_track(fig, genes, regions, 1.0, f"{contig} · {length:,} bp", length, x_range=x_range,
              color_by=color_by)
    if links is not None and len(links):
        # one tick per distinct shared segment, stacked to avoid overlap
        lk = links.assign(_r=links["discordant_rank"].where(links["is_hgt2"]))
        segs = (lk.groupby(["start", "end"])
                .agg(n=("partner", "nunique"), partners=("partner", lambda s: ", ".join(sorted(set(s))[:6])),
                     rank=("_r", lambda s: _coarsest(s)))
                .reset_index())
        by_color = {}
        for i, (_, s) in enumerate(segs.iterrows()):
            yy = 0.35 - 0.12 * (i % 3)
            color = (HGT_RAMP[_level(s["rank"])] if color_by == "hgt" and isinstance(s["rank"], str)
                     else SEGMENT)
            xs, ys = by_color.setdefault(color, ([], []))
            xs += [s["start"], s["end"], None]
            ys += [yy, yy, None]
        for color, (xs, ys) in by_color.items():
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=4),
                                     hoverinfo="skip", showlegend=False))
        fig.add_trace(go.Scatter(
            x=(segs["start"] + segs["end"]) / 2,
            y=[0.35 - 0.12 * (i % 3) for i in range(len(segs))], mode="markers",
            marker=dict(size=12, opacity=0),
            text=[f"shared with {n} contig(s)<br>{p}"
                  + (f"<br>HGT2, coarsest at {_rname(r)}" if isinstance(r, str) else "<br>same lineage (not HGT2)")
                  for n, p, r in zip(segs["n"], segs["partners"], segs["rank"])],
            hovertemplate="%{text}<extra></extra>", showlegend=False))
        fig.add_annotation(x=0, y=0.5, text="shared segments", showarrow=False, xanchor="left",
                           font=dict(size=11, color="#52514e"))
    cats = _gene_cats(genes, color_by) or {"other gene"}
    if len(regions):
        cats |= ({"kairos region"} if color_by == "hgt" else
                 {"HGT1 region" if pd.notna(x) else "other kairos region" for x in regions["hgt1_rank"]})
    _legend(fig, cats, clusters)
    _layout(fig, 260, x_range or [0, length])
    fig.update_yaxes(range=[0, 1.65])
    return fig


def _ribbons(fig, g1, g2, xf2, y2, clusters, max_pairs=400):
    """Connect genes of the same ORF cluster on the two tracks (clinker-style)."""
    m = g1[g1["cluster"] >= 0].merge(g2[g2["cluster"] >= 0], on="cluster", suffixes=("_1", "_2"))
    if m.empty:
        return
    m = m.head(max_pairs)
    by_color = {}
    for r in m.itertuples(index=False):
        a, b = sorted((xf2(r.start_2), xf2(r.end_2)))
        color = clusters.get(r.cluster, (None, SHARED_ORF, None))[1]
        xs, ys = by_color.setdefault(color, ([], []))
        xs += [r.start_1, r.end_1, b, a, r.start_1, None]
        ys += [1.0 - GENE_H, 1.0 - GENE_H, y2 + GENE_H, y2 + GENE_H, 1.0 - GENE_H, None]
    for color, (xs, ys) in by_color.items():
        fig.add_trace(go.Scatter(x=xs, y=ys, fill="toself", mode="lines", fillcolor=color,
                                 opacity=0.25, line=dict(width=0), hoverinfo="skip",
                                 showlegend=False))


def paired_map(c1, len1, genes1, regions1, c2, len2, genes2, regions2, link, pad=5000,
               color_by="hits"):
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
    if color_by != "cluster":   # in cluster mode the gene ribbons show the alignment
        fig.add_trace(go.Scatter(x=[s1, e1, b, a, s1], y=[1.0 - BAND_H, 1.0 - BAND_H, y2 + BAND_H,
                                                          y2 + BAND_H, 1.0 - BAND_H],
                                 fill="toself", fillcolor=LINK, line=dict(width=0), mode="lines",
                                 hoverinfo="skip", showlegend=False))
    clusters = {}
    if color_by == "cluster" and genes1 is not None and genes2 is not None:
        genes1, genes2 = genes1.copy(), genes2.copy()
        clusters = assign_clusters([genes1, genes2], paired=True)
        _ribbons(fig, genes1, genes2, xf2, y2, clusters)
    add_track(fig, genes1, regions1, 1.0, f"{c1} · {len1:,} bp", len1, x_range=x_range,
              color_by=color_by)
    strand_note = " (reverse complement)" if int(link.get("strand", 1) or 1) < 0 else ""
    add_track(fig, genes2, regions2, y2, f"{c2} · {len2:,} bp{strand_note}", len2, xf=xf2,
              x_range=x_range, color_by=color_by)
    cats = ({"kairos region"} if color_by == "hgt" else {"HGT1 region", "other kairos region"}) \
        | _gene_cats(genes1, color_by) \
        | _gene_cats(genes2, color_by)
    _legend(fig, cats, clusters)
    _layout(fig, 340, x_range, xtitle=f"Position on {c1} (bp)")
    fig.update_yaxes(range=[y2 - 0.55, 1.65])
    return fig
