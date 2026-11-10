import re

import pandas as pd
import streamlit as st

from hgtviz import RANK_NAMES, RANKS
from hgtviz import genemap as G

store = st.session_state["store"]
filtered: pd.DataFrame = st.session_state["filtered"]
spec = st.session_state["spec"]
rank = spec.tax_rank


def rank_name(r):
    return RANK_NAMES.get(r, "") if isinstance(r, str) else ""


st.title("Contig browser")
st.caption("Every contig with an HGT1 region or an HGT2 link. Narrow the list with the "
           "sidebar filters, then select a row to see its regions, partners and genes.")

if filtered.empty:
    st.warning("No contigs match the current filters.")
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
    "groups": filtered["groups"],
    "donor/recipient": filtered["dr_role"],
    "CARD genes": filtered["card_genes"],
    "drug classes": filtered["drug_classes"],
    "RIP-DB families": filtered["rip_families"],
})
if "insitu_as_recipient" in filtered.columns:
    tbl["in situ (rec/don)"] = (filtered["insitu_as_recipient"].astype(str) + " / "
                                + filtered["insitu_as_donor"].astype(str))
# natural sort: contig_2 before contig_10
tbl = tbl.iloc[sorted(range(len(tbl)), key=lambda i: [
    int(t) if t.isdigit() else t for t in re.split(r"(\d+)", tbl["contig"].iat[i])])]
tbl = tbl.reset_index(drop=True)

event = st.dataframe(
    tbl, hide_index=True, height=380, on_select="rerun", selection_mode="single-row",
    key="contig_table",
    column_config={
        "length": st.column_config.NumberColumn(format="localized"),
        "HGT1 regions": st.column_config.NumberColumn(width="small"),
        "HGT2 partners": st.column_config.NumberColumn(width="small"),
    },
)

# exports for the filtered set
names = set(filtered["contig"])
regions = store.regions[store.regions["contig"].isin(names)]
links = store.link_ends[store.link_ends["contig"].isin(names)]
ann = store.annotations[store.annotations["contig"].isin(names)] if store.has("annotations") else None
ex = st.columns(4)
ex[0].download_button("Contigs (TSV)", filtered.to_csv(sep="\t", index=False),
                      "contigs.tsv", "text/tab-separated-values", width="stretch")
ex[1].download_button("Regions (TSV)", regions.to_csv(sep="\t", index=False),
                      "regions.tsv", "text/tab-separated-values", width="stretch")
ex[2].download_button("Links (TSV)", links.to_csv(sep="\t", index=False),
                      "links.tsv", "text/tab-separated-values", width="stretch")
if ann is not None:
    ex[3].download_button("Annotations (TSV)", ann.to_csv(sep="\t", index=False),
                          "annotations.tsv", "text/tab-separated-values", width="stretch")

rows = event.selection.rows if event else []
if rows:
    st.session_state["sel_contig"] = tbl.iloc[rows[0]]["contig"]
sel = st.session_state.get("sel_contig")
if sel not in names:
    st.info("Select a row to see the contig's details.")
    st.stop()

# ── detail ───────────────────────────────────────────────────────────────────
c = store.contigs.set_index("contig").loc[sel]
length = int(c["length"])
st.divider()
st.subheader(sel)
lineage = " › ".join(f"{c[r]}" for r in RANKS if pd.notna(c[r]))
agree = f" · {c['tax_agreement']:.0f}% of ORFs agree" if pd.notna(c["tax_agreement"]) else ""
st.caption(f"{length:,} bp · {lineage or 'unclassified'}{agree}")

c_regions = store.regions[store.regions["contig"] == sel].sort_values("start")
c_links = store.link_ends[store.link_ends["contig"] == sel]
show_all = st.toggle("Include links between same-lineage contigs (not HGT2)", value=False,
                     help="kairos-lr also finds shared regions between contigs of the same lineage; "
                          "these are not HGT2 calls.")
if not show_all:
    c_links = c_links[c_links["is_hgt2"]]
c_dr = store.dr[(store.dr["recipient"] == sel) | (store.dr["donor"] == sel)]
c_ann = store.annotations[store.annotations["contig"] == sel] if store.has("annotations") else None
c_genes = store.genes[store.genes["contig"] == sel] if store.has("genes") else None
if c_genes is not None:
    c_genes = G.gene_category(c_genes, c_ann)

m = st.columns(6)
m[0].metric("Kairos regions", len(c_regions))
m[1].metric("HGT1 regions", int(c_regions["hgt1_rank"].notna().sum()))
m[2].metric("Partner contigs", c_links["partner"].nunique())
m[3].metric("Groups", c_regions["group_id"].nunique())
m[4].metric("CARD hits", 0 if c_ann is None else int((c_ann["source"] == "card").sum()))
m[5].metric("RIP-DB hits", 0 if c_ann is None else int((c_ann["source"] == "ripdb").sum()))

# gene map with optional zoom to a region
zoom_opts = ["Whole contig"] + c_regions["region_idx"].tolist()
zoom = st.selectbox("Zoom to", zoom_opts, key=f"zoom_{sel}")
x_range = None
if zoom != "Whole contig":
    r = c_regions.set_index("region_idx").loc[zoom]
    pad = max(2000, int((r["end"] - r["start"]) * 0.3))
    x_range = [max(0, r["start"] - pad), min(length, r["end"] + pad)]
st.plotly_chart(G.contig_map(sel, length, c_genes, c_regions, c_links, x_range=x_range),
                width="stretch", key=f"map_{sel}_{zoom}")
if c_genes is None:
    st.caption("No gene coordinates in this store (ffn_coords.tsv was not found at build time).")

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
        pt = c_links.copy()
        pt["lineages differ at"] = pt["discordant_rank"].map(rank_name)
        cols = ["partner", f"p_{rank}", "lineages differ at", "region", "partner_region",
                "start", "end", "partner_start", "partner_end", "strand", "identity",
                "n_shared_orfs", "is_hgt2"]
        cols = [x for x in cols if x in pt.columns]
        st.dataframe(pt[cols].rename(columns={f"p_{rank}": f"partner {RANK_NAMES[rank].lower()}"}),
                     hide_index=True, width="stretch")
        st.markdown("**Compare with a partner**")
        labels = [f"{p} · {rg} ↔ {pr}" for p, rg, pr in
                  zip(pt["partner"], pt["region"], pt["partner_region"])]
        choice = st.selectbox("Link", ["—"] + labels, key=f"pair_{sel}", label_visibility="collapsed")
        if choice != "—":
            link = pt.iloc[labels.index(choice)]
            p = link["partner"]
            p_len = store.regions.loc[store.regions["contig"] == p, "contig_length"]
            p_len = int(p_len.iloc[0]) if len(p_len) else int(max(link["partner_end"], 1))
            p_regions = store.regions[store.regions["contig"] == p]
            p_ann = store.annotations[store.annotations["contig"] == p] if store.has("annotations") else None
            p_genes = store.genes[store.genes["contig"] == p] if store.has("genes") else None
            if p_genes is not None:
                p_genes = G.gene_category(p_genes, p_ann)
            pad = st.slider("Flank (bp)", 0, 50000, 5000, 1000, key=f"pad_{sel}")
            st.plotly_chart(G.paired_map(sel, length, c_genes, c_regions, p, p_len, p_genes,
                                         p_regions, link, pad=pad),
                            width="stretch", key=f"pairfig_{sel}_{choice}")

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
    if c_ann is None or c_ann.empty:
        st.caption("No CARD or RIP-DB hits on this contig.")
    else:
        st.dataframe(c_ann[["orf_id", "source", "label", "accession", "family", "category",
                            "mechanism", "pident", "bitscore", "evalue", "region_idx"]]
                     .sort_values(["source", "orf_id"]),
                     hide_index=True, width="stretch")

if store.has("insitu"):
    with tabs[5]:
        ins = store.insitu
        ins = ins[(ins["recipient_contig"] == sel) | (ins["donor_contig"] == sel)]
        if ins.empty:
            st.caption("No in situ calls involve this contig.")
        else:
            st.dataframe(ins, hide_index=True, width="stretch")
