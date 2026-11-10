"""hgt_viz entry point: dataset + filters in the sidebar, pages via st.navigation."""
import os

import streamlit as st

from hgtviz import RANK_NAMES, RANKS
from hgtviz import filters as F
from hgtviz import store as S
from hgtviz.build import build

st.set_page_config(page_title="HGT viz", layout="wide")


@st.cache_resource(show_spinner="Loading store…", max_entries=4)
def load_store(path: str, stamp: float) -> S.Store:
    return S.load(path)


def _stamp(path):
    p = os.path.join(path, "meta.json")
    return os.path.getmtime(p) if os.path.exists(p) else 0.0


def _looks_like_results(path):
    return os.path.isdir(os.path.join(path, "detect")) or os.path.isdir(os.path.join(path, "results", "detect"))


def dataset_panel():
    st.sidebar.header("Dataset")
    if "store_path" not in st.session_state:
        st.session_state["store_path"] = os.environ.get("HGT_VIZ_STORE", os.path.join(os.getcwd(), "store"))
    path = st.sidebar.text_input("Store directory", key="store_path",
                                 help="Built with `python -m hgtviz.build --results … --out …`")
    if S.is_store(path):
        return load_store(path, _stamp(path))

    st.sidebar.warning("No hgt_viz store here.")
    if os.path.isdir(path) and _looks_like_results(path):
        results = path if os.path.isdir(os.path.join(path, "detect")) else os.path.join(path, "results")
        out = os.path.join(path, "hgt_viz_store")
        card = st.sidebar.text_input("CARD aro_index.tsv (optional)")
        if st.sidebar.button(f"Build store into {out}"):
            msgs = []
            with st.spinner("Building…"):
                build(results, out, card_index=card or None, log=msgs.append)
            st.sidebar.code("\n".join(msgs))
            st.session_state["store_path"] = out
            st.rerun()
    return None


def _reset_taxa():
    for k in ["f_contig_taxa", "f_region_taxa", "f_partner_taxa"]:
        st.session_state[k] = []


def filter_panel(store: S.Store) -> F.FilterSpec:
    sb = st.sidebar
    sb.header("Filter contigs")
    if sb.button("Clear all filters", width="stretch"):
        for k in list(st.session_state):
            if k.startswith("f_"):
                del st.session_state[k]
        st.rerun()
    spec = F.FilterSpec()
    spec.search = sb.text_input("Contig ID contains", key="f_search")
    evidence = ["HGT1", "HGT2", "Donor–recipient"] + (["In situ"] if store.has("insitu") else [])
    spec.require = sb.multiselect("Require evidence (all of)", evidence, key="f_require")

    with sb.expander("Taxonomy", expanded=True):
        spec.tax_rank = st.selectbox("Rank", RANKS, index=1, format_func=RANK_NAMES.get,
                                     key="f_tax_rank", on_change=_reset_taxa)
        opts = F.options(store, spec.tax_rank)
        spec.contig_taxa = st.multiselect("Host contig taxon", opts["contig_taxa"], key="f_contig_taxa")
        spec.region_taxa = st.multiselect("HGT1 region taxon (putative donor lineage)",
                                          opts["region_taxa"], key="f_region_taxa")
        spec.partner_taxa = st.multiselect("HGT2 partner contig taxon", opts["partner_taxa"],
                                           key="f_partner_taxa")

    with sb.expander("HGT calls"):
        spec.hgt1_ranks = st.multiselect("HGT1 detected at", RANKS, format_func=RANK_NAMES.get,
                                         key="f_hgt1_ranks")
        spec.hgt2_ranks = st.multiselect("HGT2 lineages differ at", RANKS, format_func=RANK_NAMES.get,
                                         key="f_hgt2_ranks")
        spec.min_partners = st.number_input("Min. HGT2 partner contigs", 0, step=1, key="f_min_partners")
        spec.groups = st.multiselect("HGT group", opts["groups"], key="f_groups")
        spec.min_group_size = st.number_input("Min. group size (regions)", 0, step=1, key="f_min_group")
        spec.min_length = st.number_input("Min. contig length (bp)", 0, step=1000, key="f_min_length")

    with sb.expander("Function", expanded=bool(opts["card_genes"] or opts["rip_families"])):
        if store.annotations is None:
            st.caption("No CARD / RIP-DB annotations in this store.")
        scope = st.radio("Where", ["anywhere on the contig", "inside a kairos region"],
                         key="f_scope", horizontal=False)
        spec.function_scope = "region" if scope.startswith("inside") else "contig"
        spec.card_genes = st.multiselect("CARD gene", opts["card_genes"], key="f_card")
        spec.drug_classes = st.multiselect("CARD drug class", opts["drug_classes"], key="f_drug",
                                           help="Needs --card-index at build time")
        spec.rip_families = st.multiselect("RIP-DB family", opts["rip_families"], key="f_rip")
    return spec


store = dataset_panel()
pages = [
    st.Page("views/browser.py", title="Contig browser", icon=":material/table_rows:", default=True),
    st.Page("views/summary.py", title="Summary", icon=":material/bar_chart:"),
    st.Page("views/groups.py", title="HGT groups", icon=":material/hub:"),
]
pg = st.navigation(pages)

if store is None:
    st.title("HGT viz")
    st.info("Point the sidebar at a store directory built from kairos-lr results:\n\n"
            "```\npython -m hgtviz.build --results path/to/results --out path/to/store "
            "--card-index path/to/aro_index.tsv\n```")
    st.stop()

spec = filter_panel(store)
filtered = F.apply(store, spec)
st.session_state["store"] = store
st.session_state["spec"] = spec
st.session_state["filtered"] = filtered
st.sidebar.caption(f"**{len(filtered):,}** of {len(store.contigs):,} contigs · "
                   f"catalogue `{store.meta.get('catalogue')}` · built {store.meta.get('built', '?')}")
pg.run()
