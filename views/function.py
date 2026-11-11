"""Search genes by function; the matching contigs become a sidebar filter."""
import pandas as pd
import streamlit as st

from hgtviz import HGT1_NAME, HGT2_NAME, RANK_NAMES, RANKS
from hgtviz import cog as COG
from hgtviz import filters as F

store = st.session_state["store"]
filtered: pd.DataFrame = st.session_state["filtered"]
spec: F.FilterSpec = st.session_state["spec"]
LIMIT = 5000

st.title("Function search")
st.caption("Find genes by annotation. **Apply** turns the search into a contig filter: only "
           "contigs that carry a match stay in the browser, summary and groups pages. Each "
           "criterion you set must be met somewhere on the contig; the gene table lists every "
           "gene that meets any of them.")

opts = F.options(store, spec.tax_rank)
cur = st.session_state.get("f_func") or {}
has_ann = store.has("annotations")
has_egg = store.has("eggnog")
if not (has_ann or has_egg):
    st.warning("This store has no CARD / RIP-DB or eggNOG annotations.")
    st.stop()

with st.form("function_form", border=True):
    c = st.columns(3)
    card = c[0].multiselect("CARD gene", opts["card_genes"], default=cur.get("card_genes", []),
                            disabled=not has_ann)
    drug = c[1].multiselect("CARD drug class", opts["drug_classes"],
                            default=cur.get("drug_classes", []), disabled=not has_ann,
                            help="Needs --card-index at build time")
    rip = c[2].multiselect("RIP-DB family", opts["rip_families"], default=cur.get("rip_families", []),
                           disabled=not has_ann)
    c = st.columns([2, 1])
    egg = c[0].text_input("eggNOG annotation contains", value=cur.get("eggnog_text", ""),
                          disabled=not has_egg, placeholder="e.g. transposase, K07483, Rep_3, recombinase",
                          help="Case-insensitive match on the eggNOG preferred name, description, "
                               "KEGG KO, Pfam, CAZy and EC of a gene.")
    cogs = c[1].multiselect("COG category (eggNOG)", opts["cog_letters"], format_func=COG.label,
                            default=cur.get("cog_letters", []), disabled=not has_egg)
    scopes = F.SCOPES
    c = st.columns([2, 1])
    scope = c[0].radio("Where", list(scopes), format_func=scopes.get, horizontal=True,
                       index=list(scopes).index(cur.get("function_scope", "contig")),
                       help=f"'Region with an HGT call': the gene lies in a kairos region that has a "
                            f"{HGT1_NAME.lower()} (HGT1) call or a {HGT2_NAME.lower()} (HGT2) link.")
    hgt_ranks = c[1].multiselect("…with HGT at (coarsest rank)", RANKS, format_func=RANK_NAMES.get,
                                 default=cur.get("function_hgt_ranks", []),
                                 help="Only for 'inside a region with an HGT call'. Empty = any rank. "
                                      "The rank is the coarsest over the region's HGT1 call and HGT2 links.")
    b = st.columns([1, 1, 4])
    apply = b[0].form_submit_button("Apply as filter", type="primary", width="stretch")
    clear = b[1].form_submit_button("Clear", width="stretch")

if apply:
    st.session_state["f_func"] = {"card_genes": card, "drug_classes": drug, "rip_families": rip,
                                  "cog_letters": cogs, "eggnog_text": egg.strip(),
                                  "function_scope": scope,
                                  "function_hgt_ranks": hgt_ranks if scope == "hgt" else []}
    st.rerun()
if clear:
    st.session_state.pop("f_func", None)
    st.rerun()

if not F.function_active(spec):
    st.info("Set at least one criterion and press **Apply as filter**.")
    st.stop()

# ── results ──────────────────────────────────────────────────────────────────
with st.spinner("Searching annotations…"):
    all_hits, all_tot = F.function_hits(store, spec, None, LIMIT)
    hits, tot = F.function_hits(store, spec, set(filtered["contig"]), LIMIT)

m = st.columns(3)
m[0].metric("Matching genes (whole store)", f"{sum(all_tot.values()):,}",
            help=" · ".join(f"{k}: {v:,}" for k, v in all_tot.items()))
m[1].metric("Matching genes on filtered contigs", f"{sum(tot.values()):,}",
            help=" · ".join(f"{k}: {v:,}" for k, v in tot.items()))
m[2].metric("Contigs passing all filters", f"{len(filtered):,}")

if not sum(all_tot.values()):
    st.warning("No gene in this store matches. For eggNOG text, try a shorter term "
               "(e.g. 'transpos' rather than 'transposase IS200').")
    st.stop()
if filtered.empty:
    st.warning("Genes match, but no contig meets every criterion together with the other sidebar "
               "filters. Criteria are combined with AND: a contig needs, for example, both the CARD "
               "gene and the eggNOG text. Clear some sidebar filters or criteria.")
    st.stop()
if any(v > LIMIT for v in tot.values()):
    st.caption(f"Showing at most {LIMIT:,} genes per source.")

def rank_name(r):
    return RANK_NAMES.get(r, "") if isinstance(r, str) else ""


in_hgt = hits["hgt_rank"].notna()
st.markdown(f"Of the matching genes on the filtered contigs, **{int(in_hgt.sum()):,}** lie in a region "
            f"with an HGT call, **{int((hits['in_region'] & ~in_hgt).sum()):,}** in a kairos region "
            f"without one and **{int((~hits['in_region'].astype(bool)).sum()):,}** outside kairos regions"
            + (" (in the table below)." if any(v > LIMIT for v in tot.values()) else "."))
only_hgt = st.toggle("Only genes in regions with an HGT call", value=False, key="fn_only_hgt")
if only_hgt:
    hits = hits[in_hgt]

c = store.contigs.set_index("contig")
r = spec.tax_rank
tbl = hits.assign(
    **{"HGT level": hits["hgt_rank"].map(rank_name),
       "HGT1 at": hits["hgt1_rank"].map(rank_name),
       "HGT2 at": hits["hgt2_rank"].map(rank_name),
       "HGT2 partners": hits["hgt2_partners"],
       RANK_NAMES[r]: hits["contig"].map(c[r])})
tbl = tbl.drop(columns=["hgt_rank", "hgt1_rank", "hgt2_rank", "hgt2_partners"])
tbl = tbl.rename(columns={"orf_id": "gene", "region_idx": "kairos region", "in_region": "in region",
                          "source": "matched", "detail": "description"})
order = ["contig", "gene", "matched", "name", "description", "kairos region", "HGT level",
         "HGT1 at", "HGT2 at", "HGT2 partners", RANK_NAMES[r], "in region"]
tbl = tbl[order].sort_values(["contig", "gene"], key=lambda s: s.astype(str)).reset_index(drop=True)
st.caption("Select a gene to open its contig in the browser.")
ev = st.dataframe(tbl, hide_index=True, width="stretch", height=420, on_select="rerun",
                  selection_mode="single-row", key="function_hits",
                  column_config={
                      "description": st.column_config.TextColumn(width="large"),
                      "HGT level": st.column_config.TextColumn(
                          help="Coarsest rank of the HGT evidence on the gene's region (HGT1 or HGT2); "
                               "blank = no HGT call on the region"),
                      "HGT1 at": st.column_config.TextColumn(
                          help=f"{HGT1_NAME}: rank at which the region's taxonomy differs from its contig"),
                      "HGT2 at": st.column_config.TextColumn(
                          help=f"{HGT2_NAME}: coarsest rank at which the region's partner contigs differ"),
                      "HGT2 partners": st.column_config.NumberColumn(
                          help="Contigs of another lineage that share this region")})
rows = ev.selection.rows if ev else []
if rows:
    g = tbl.iloc[rows[0]]
    if st.button(f"Open {g['contig']} in the contig browser", type="primary"):
        st.session_state["sel_contig"] = g["contig"]
        st.switch_page("views/browser.py")
