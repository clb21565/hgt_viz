"""Search genes by function; the matching contigs become a sidebar filter."""
import pandas as pd
import streamlit as st

from hgtviz import HGT1_NAME, HGT2_NAME, RANK_NAMES
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
    scopes = {"contig": "anywhere on the contig", "region": "inside a kairos region"}
    scope = st.radio("Where", list(scopes), format_func=scopes.get, horizontal=True,
                     index=list(scopes).index(cur.get("function_scope", "contig")))
    b = st.columns([1, 1, 4])
    apply = b[0].form_submit_button("Apply as filter", type="primary", width="stretch")
    clear = b[1].form_submit_button("Clear", width="stretch")

if apply:
    st.session_state["f_func"] = {"card_genes": card, "drug_classes": drug, "rip_families": rip,
                                  "cog_letters": cogs, "eggnog_text": egg.strip(),
                                  "function_scope": scope}
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

c = store.contigs.set_index("contig")
r = spec.tax_rank
tbl = hits.assign(
    **{RANK_NAMES[r]: hits["contig"].map(c[r]),
       f"{HGT1_NAME} regions": hits["contig"].map(c["n_hgt1_regions"]),
       f"{HGT2_NAME} partners": hits["contig"].map(c["n_hgt2_partners"])})
tbl = tbl.rename(columns={"orf_id": "gene", "region_idx": "kairos region", "in_region": "in region",
                          "source": "matched", "detail": "description"})
tbl = tbl.sort_values(["contig", "gene"], key=lambda s: s.astype(str)).reset_index(drop=True)
st.caption("Select a gene to open its contig in the browser.")
ev = st.dataframe(tbl, hide_index=True, width="stretch", height=420, on_select="rerun",
                  selection_mode="single-row", key="function_hits",
                  column_config={"description": st.column_config.TextColumn(width="large")})
rows = ev.selection.rows if ev else []
if rows:
    g = tbl.iloc[rows[0]]
    if st.button(f"Open {g['contig']} in the contig browser", type="primary"):
        st.session_state["sel_contig"] = g["contig"]
        st.switch_page("views/browser.py")
