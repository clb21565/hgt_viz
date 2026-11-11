import pandas as pd
import streamlit as st

from hgtviz import HGT1_LABEL, HGT1_NAME, HGT2_LABEL, HGT2_NAME, RANK_NAMES, RANKS
from hgtviz import genemap as G
from hgtviz import lineages as LIN

store = st.session_state["store"]
filtered: pd.DataFrame = st.session_state["filtered"]
spec = st.session_state["spec"]
rank = spec.tax_rank


def rank_name(r):
    return RANK_NAMES.get(r, "") if isinstance(r, str) else ""


st.title("Contig browser")
st.caption(f"Every contig with a {HGT1_NAME.lower()} (HGT1) or a {HGT2_NAME.lower()} (HGT2). "
           "Narrow the list with the sidebar filters, then select a row to see its regions, "
           "partners and genes. The 'HGT calls explained' page describes both types.")

if filtered.empty:
    st.warning("No contigs match the current filters.")
    st.stop()

QUICK = {"any": "Any HGT", "hgt1": HGT1_LABEL, "hgt2": HGT2_LABEL, "both": "Both",
         "hgt1_only": f"{HGT1_NAME} only", "hgt2_only": f"{HGT2_NAME} only"}
quick = st.segmented_control("Show contigs with", list(QUICK), format_func=QUICK.get,
                             default="any", key="quick_evidence") or "any"
h1, h2 = filtered["hgt1"], filtered["hgt2"]
filtered = filtered[{"any": h1 | h2, "hgt1": h1, "hgt2": h2, "both": h1 & h2,
                     "hgt1_only": h1 & ~h2, "hgt2_only": h2 & ~h1}[quick]]
if filtered.empty:
    st.warning("No contigs of this kind pass the current filters.")
    st.stop()

# ── table ────────────────────────────────────────────────────────────────────
tbl = pd.DataFrame({
    "contig": filtered["contig"],
    "length": filtered["length"],
    RANK_NAMES[rank]: filtered[rank],
    "taxon": filtered["taxon"],
    "HGT1 regions": filtered["n_hgt1_regions"],
    "HGT1 rank": filtered["hgt1_rank"].map(rank_name),
    "HGT2 partners": filtered["n_hgt2_partners"],
    "HGT2 rank": filtered["hgt2_rank"].map(rank_name),
    "lineage": filtered[RANKS].astype("object").where(filtered[RANKS].notna(), "")
               .astype(str).agg("; ".join, axis=1).str.replace(r"(; )+$", "", regex=True),
    "groups": filtered["groups"],
    "donor/recipient": filtered["dr_role"],
    "CARD genes": filtered["card_genes"],
    "drug classes": filtered["drug_classes"],
    "RIP-DB families": filtered["rip_families"],
})
if "insitu_as_recipient" in filtered.columns:
    tbl["in situ (rec/don)"] = (filtered["insitu_as_recipient"].astype(str) + " / "
                                + filtered["insitu_as_donor"].astype(str))
tbl = tbl.reset_index(drop=True)  # store.contigs is natural-sorted at load

event = st.dataframe(
    tbl, hide_index=True, height=380, on_select="rerun", selection_mode="single-row",
    key="contig_table",
    column_config={
        "length": st.column_config.NumberColumn(format="localized"),
        "HGT1 regions": st.column_config.NumberColumn(
            width="small", help=f"{HGT1_NAME} regions: region taxonomy differs from this contig's"),
        "HGT1 rank": st.column_config.TextColumn(help=f"{HGT1_NAME}: coarsest rank of disagreement"),
        "HGT2 partners": st.column_config.NumberColumn(
            width="small", help=f"{HGT2_NAME}: contigs of another lineage sharing a region with this one"),
        "HGT2 rank": st.column_config.TextColumn(help=f"{HGT2_NAME}: coarsest rank at which the lineages differ"),
        "lineage": st.column_config.TextColumn(width="medium"),
    },
)

# exports for the filtered set: the TSV is only written when a button is clicked
names = set(filtered["contig"])


def _tsv(get):
    return lambda: get().to_csv(sep="\t", index=False)


exports = [
    ("Contigs", lambda: filtered),
    ("Regions", lambda: store.regions[store.regions["contig"].isin(names)]),
    ("Links", lambda: store.link_ends[store.link_ends["contig"].isin(names)]),
]
if store.has("annotations"):
    exports.append(("Annotations", lambda: store.annotations[store.annotations["contig"].isin(names)]))
ex = st.columns(4)
for col, (label, get) in zip(ex, exports):
    col.download_button(f"{label} (TSV)", _tsv(get), f"{label.lower()}.tsv",
                        "text/tab-separated-values", width="stretch", on_click="ignore")

rows = event.selection.rows if event else []
if rows:
    st.session_state["sel_contig"] = tbl.iloc[rows[0]]["contig"]
sel = st.session_state.get("sel_contig")
if sel not in names:
    st.info("Select a row to see the contig's details.")
    st.stop()

# ── detail ───────────────────────────────────────────────────────────────────
# A fragment: zoom, tabs, partner choice and sliders rerun only this part, not the
# sidebar filters and the contig table.
@st.fragment
def detail(sel):
    _detail(sel)


def _detail(sel):
    c = store.rows_for("contigs", sel).iloc[0]
    length = int(c["length"])
    st.divider()
    st.subheader(sel)
    lineage = " › ".join(f"{c[r]}" for r in RANKS if pd.notna(c[r]))
    agree = f" · {c['tax_agreement']:.0f}% of ORFs agree" if pd.notna(c["tax_agreement"]) else ""
    st.caption(f"{length:,} bp · {lineage or 'unclassified'}{agree}")

    c_regions = store.rows_for("regions", sel).sort_values("start")
    c_links = store.rows_for("link_ends", sel)
    show_all = st.toggle(f"Include links between same-lineage contigs (not {HGT2_NAME.lower()})",
                         value=False,
                         help="kairos-lr also finds shared regions between contigs of the same lineage; "
                              "these are not HGT2 calls.")
    if not show_all:
        c_links = c_links[c_links["is_hgt2"]]
    c_dr = store.dr[(store.dr["recipient"] == sel) | (store.dr["donor"] == sel)]
    c_ann = store.rows_for("annotations", sel) if store.has("annotations") else None
    c_genes = store.rows_for("genes", sel) if store.has("genes") else None
    c_egg = store.eggnog_for(sel) if store.has("eggnog") else None
    if c_genes is not None:
        c_genes = G.gene_category(c_genes, c_ann, c_egg)
    c_links_all = store.rows_for("link_ends", sel)

    m = st.columns(6)
    m[0].metric("Kairos regions", len(c_regions))
    m[1].metric(f"{HGT1_NAME} regions", int(c_regions["hgt1_rank"].notna().sum()))
    m[2].metric("Partner contigs", c_links["partner"].nunique())
    m[3].metric("Groups", c_regions["group_id"].nunique())
    m[4].metric("CARD hits", 0 if c_ann is None else int((c_ann["source"] == "card").sum()))
    m[5].metric("RIP-DB hits", 0 if c_ann is None else int((c_ann["source"] == "ripdb").sum()))

    # gene map with optional zoom to a region
    zoom_opts = ["Whole contig"] + c_regions["region_idx"].tolist()
    pending = st.session_state.pop("zoom_pending", None)   # set by the Function search page
    if pending in zoom_opts:
        st.session_state[f"zoom_{sel}"] = pending
    hl = st.session_state.get("sel_gene")
    if c_genes is None or hl not in set(c_genes["orf_id"]):
        hl = None
    if hl:
        g = c_genes[c_genes["orf_id"] == hl].iloc[0]
        nm = f" ({g['gene_name']})" if g.get("gene_name") else ""
        hc = st.columns([5, 1])
        hc[0].markdown(f":red[**▬**] Selected gene **{hl}**{nm}, {int(g['start']) + 1:,}–{int(g['end']):,} bp"
                       + (f", in {g['region_idx']}" if isinstance(g.get("region_idx"), str) else ""))
        hc[1].button("Clear", key="clear_sel_gene", width="stretch",
                     on_click=lambda: st.session_state.pop("sel_gene", None))
    zc = st.columns([2, 3])
    zoom = zc[0].selectbox("Zoom to", zoom_opts, key=f"zoom_{sel}")
    color_opts = {"HGT evidence": "hgt", "CARD / RIP-DB hits": "hits"}
    if store.has("eggnog"):
        color_opts["COG group (eggNOG)"] = "cog"
    if store.has("genes") and "cluster" in store.genes.columns:
        color_opts["Identical ORFs (mmseqs clusters)"] = "cluster"
    color_by = color_opts[zc[1].radio("Colour genes by", list(color_opts), horizontal=True,
                                      key="color_by")]
    oc = st.columns([3, 1])
    EV = {"both": "Both", "hgt1": HGT1_LABEL, "hgt2": HGT2_LABEL}
    ev_sel = oc[0].segmented_control(
        "Evidence shown", list(EV), format_func=EV.get, default="both", key="evidence_shown",
        help="Limit the HGT colouring, region labels and shared-segment ticks to one call type.") or "both"
    use = {"both": ("hgt1", "hgt2"), "hgt1": ("hgt1",), "hgt2": ("hgt2",)}[ev_sel]
    names_on = oc[1].toggle("Gene names", value=True, key="gene_names",
                            help="CARD gene or eggNOG preferred name under each gene, where it fits. "
                                 "Zoom in (Zoom to) to fit more.")
    c_genes, c_regions_h = G.annotate_hgt(c_genes, c_regions, c_links_all, c, use=use)
    x_range = None
    if zoom != "Whole contig":
        r = c_regions.set_index("region_idx").loc[zoom]
        pad = max(2000, int((r["end"] - r["start"]) * 0.3))
        x_range = [max(0, r["start"] - pad), min(length, r["end"] + pad)]
    st.plotly_chart(G.contig_map(sel, length, c_genes, c_regions_h,
                                 None if ev_sel == "hgt1" else c_links, x_range=x_range,
                                 color_by=color_by, names=names_on, highlight=hl),
                    width="stretch", key=f"map_{sel}_{zoom}_{color_by}_{ev_sel}_{names_on}_{hl}")
    if c_genes is None:
        st.caption("No gene coordinates in this store (ffn_coords.tsv was not found at build time).")

    # ── lineages involved ────────────────────────────────────────────────────
    lc = st.columns([3, 1])
    full = lc[1].toggle("Full partner lineages", value=False, key="lin_full",
                        help="Off: partners are grouped by the taxon at which they differ from this "
                             "contig (lower ranks hidden). On: one row per distinct full lineage.")
    lin, hidden = LIN.hgt_lineages(c, c_regions, c_links_all, use=use, full=full)
    if len(lin) > 1:
        lc[0].markdown("**Lineages in this contig's HGT calls**")
        st.caption("Host first, then each foreign region's own classification and each distinct "
                   "partner lineage. Grey = same as the host; dark orange = the rank where it first "
                   "differs; light orange = also differs below that."
                   + (f" {hidden:,} rarer partner lineages not shown (Partners tab has all)." if hidden else ""))
        st.dataframe(LIN.style(lin), hide_index=True, width="stretch",
                     height=min(36 * (len(lin) + 1) + 3, 420),
                     column_config={"contigs": st.column_config.TextColumn(
                         help="Partner contigs with this lineage", alignment="right")})

    # ── tabs ─────────────────────────────────────────────────────────────────────
    tabs = st.tabs(["Regions", f"Partners ({c_links['partner'].nunique()})",
                    f"Donor–recipient ({len(c_dr)})", "Groups", "Annotations"]
                   + (["In situ"] if store.has("insitu") else []))

    with tabs[0]:
        cols = ["region_idx", "start", "end", "region_len", "region_type", "hgt1_rank",
                "region_taxon"] + [f"r_{r}" for r in RANKS] + ["group_id", "group_size"]
        rt = c_regions[cols].assign(hgt1_rank=c_regions["hgt1_rank"].map(rank_name))
        st.dataframe(rt.rename(columns={f"r_{r}": f"region {RANK_NAMES[r].lower()}"
                                                     for r in RANKS}),
                     hide_index=True, width="stretch")

    with tabs[1]:
        if c_links.empty:
            st.caption("No partner contigs.")
        else:
            pt = c_links.reset_index(drop=True)
            pt["lineages differ at"] = pt["discordant_rank"].map(rank_name)
            cols = ["partner", f"p_{rank}", "lineages differ at", "region", "partner_region",
                    "start", "end", "partner_start", "partner_end", "strand", "identity",
                    "n_shared_orfs", "is_hgt2"]
            cols = [x for x in cols if x in pt.columns]
            st.caption("Select a row to compare the two contigs around the shared segment.")
            ev = st.dataframe(pt[cols].rename(columns={f"p_{rank}": f"partner {RANK_NAMES[rank].lower()}"}),
                              hide_index=True, width="stretch", height=300, on_select="rerun",
                              selection_mode="single-row", key=f"pairs_{sel}_{show_all}")
            prow = ev.selection.rows if ev else []
            if prow:
                link = pt.iloc[prow[0]]
                choice = f"{link['partner']}_{link['region']}_{link['partner_region']}"
                p = link["partner"]
                p_len = store.rows_for("regions", p)["contig_length"]
                p_len = int(p_len.iloc[0]) if len(p_len) else int(max(link["partner_end"], 1))
                p_regions = store.rows_for("regions", p)
                p_ann = store.rows_for("annotations", p) if store.has("annotations") else None
                p_genes = store.rows_for("genes", p) if store.has("genes") else None
                if p_genes is not None:
                    p_egg = store.eggnog_for(p) if store.has("eggnog") else None
                    p_genes = G.gene_category(p_genes, p_ann, p_egg)
                p_host = store.rows_for("contigs", p)
                p_genes, p_regions = G.annotate_hgt(p_genes, p_regions, store.rows_for("link_ends", p),
                                                    p_host.iloc[0] if len(p_host) else None, use=use)
                pad = st.slider("Flank (bp)", 0, 50000, 5000, 1000, key=f"pad_{sel}")
                st.plotly_chart(G.paired_map(sel, length, c_genes, c_regions_h, p, p_len, p_genes,
                                             p_regions, link, pad=pad, color_by=color_by,
                                             names=names_on, highlight=hl),
                                width="stretch",
                                key=f"pairfig_{sel}_{choice}_{color_by}_{ev_sel}_{names_on}_{hl}")

    with tabs[2]:
        if c_dr.empty:
            st.caption("Not part of any donor–recipient pair.")
        else:
            d = c_dr.copy()
            d["role"] = d["recipient"].eq(sel).map({True: "recipient", False: "donor"})
            d["rank"] = d["rank"].map(rank_name)
            st.dataframe(d[["role", "recipient", "recipient_region", "donor", "donor_region", "rank"]],
                         hide_index=True, width="stretch")

    with tabs[3]:
        gids = c_regions["group_id"].dropna().astype(int).unique()
        if len(gids) == 0:
            st.caption("No regions of this contig are in an HGT group.")
        else:
            if store.groups is not None:
                st.dataframe(store.groups[store.groups["hgt_group_id"].isin(gids)],
                             hide_index=True, width="stretch")
            members = store.regions[store.regions["group_id"].isin(gids)]
            st.caption(f"{len(members):,} member regions on {members['contig'].nunique():,} contigs. "
                       "Open the HGT groups page for the full member list.")

    with tabs[4]:
        st.markdown("**CARD / RIP-DB**")
        if c_ann is None or c_ann.empty:
            st.caption("No CARD or RIP-DB hits on this contig.")
        else:
            st.dataframe(c_ann[["orf_id", "source", "label", "accession", "family", "category",
                                "mechanism", "pident", "bitscore", "evalue", "region_idx"]]
                         .sort_values(["source", "orf_id"]),
                         hide_index=True, width="stretch",
                         column_config={"evalue": st.column_config.NumberColumn(format="%.1e")})
        if store.has("eggnog"):
            st.markdown("**eggNOG**")
            if c_egg is None or c_egg.empty:
                st.caption("No eggNOG annotations on this contig.")
            else:
                only = st.toggle("Only genes inside kairos regions", value=True, key=f"egg_reg_{sel}")
                e = c_egg[c_egg["in_region"]] if only else c_egg
                e = e.assign(_n=pd.to_numeric(e["orf_id"].str.rsplit("_", n=1).str[-1], errors="coerce")).sort_values("_n")
                tcols = ["orf_id", "region_idx", "preferred_name", "description", "cog_category",
                         "kegg_ko", "pfams", "cazy", "ec"]
                e = e[tcols + ["evalue"]].astype({c: "object" for c in tcols})
                e[tcols] = e[tcols].where(e[tcols].notna(), "")
                st.dataframe(e,
                             hide_index=True, width="stretch",
                             column_config={"evalue": st.column_config.NumberColumn(format="%.1e")})
                st.caption(f"{len(c_egg):,} of {len(c_genes) if c_genes is not None else '?'} "
                           f"genes on this contig have an eggNOG hit.")

    if store.has("insitu"):
        with tabs[5]:
            ins = store.insitu
            ins = ins[(ins["recipient_contig"] == sel) | (ins["donor_contig"] == sel)]
            if ins.empty:
                st.caption("No in situ calls involve this contig.")
            else:
                st.dataframe(ins, hide_index=True, width="stretch")


detail(sel)
