import pandas as pd
import streamlit as st

from hgtviz import RANK_NAMES, RANKS

store = st.session_state["store"]
filtered: pd.DataFrame = st.session_state["filtered"]

st.title("HGT groups")
st.caption("Groups from get_hgt-v2.py (union-find over HGT2 region pairs) that contain at least "
           "one region on a contig in the filtered set.")

R = store.regions
if R["group_id"].isna().all():
    st.info("This store has no groups (no hgt_groups.tsv at build time).")
    st.stop()

names = set(filtered["contig"])
touch = store.subset("regions", filtered)
touch = touch[touch["group_id"].notna()]
if touch.empty:
    st.warning("No groups touch the filtered contigs.")
    st.stop()

gids = touch["group_id"].astype(int).unique()
if store.groups is not None:
    g = store.groups[store.groups["hgt_group_id"].isin(gids)].copy()
else:
    g = (R[R["group_id"].isin(gids)].groupby("group_id")
         .agg(num_regions=("region_idx", "nunique"), num_contigs=("contig", "nunique"))
         .reset_index().rename(columns={"group_id": "hgt_group_id"}))
in_set = touch.groupby("group_id")["contig"].nunique()
g.insert(1, "contigs_in_filter", g["hgt_group_id"].map(in_set).fillna(0).astype(int))

# function content of each group's member regions
if store.has("annotations"):
    A = store.annotations[store.annotations["in_region"]]
    A = A.merge(R[["region_idx", "group_id"]], on="region_idx")
    for src, col in [("card", "CARD genes"), ("ripdb", "RIP-DB families")]:
        s = A[A["source"] == src].groupby("group_id")["label"].agg(
            lambda x: ", ".join(f"{k} ({v})" for k, v in x.value_counts().items()))
        g[col] = g["hgt_group_id"].map(s).fillna("")

g = g.sort_values("contigs_in_filter", ascending=False).reset_index(drop=True)
event = st.dataframe(g, hide_index=True, height=360, on_select="rerun",
                     selection_mode="single-row", key="group_table")
rows = event.selection.rows if event else []
if not rows:
    st.info("Select a group to list its members.")
    st.stop()

gid = int(g.iloc[rows[0]]["hgt_group_id"])
members = R[R["group_id"] == gid].merge(
    store.contigs[["contig"] + RANKS], on="contig", how="left", suffixes=("", "_contig"))
st.subheader(f"Group {gid}")
st.caption(f"{len(members):,} regions on {members['contig'].nunique():,} contigs. "
           "Card/RIP columns count hits inside each member region; contig ranks are the host contig's lineage.")
if store.has("annotations"):
    A = store.annotations[store.annotations["region_idx"].isin(members["region_idx"])]
    for src, col in [("card", "CARD"), ("ripdb", "RIP-DB")]:
        members[col] = members["region_idx"].map(
            A[A["source"] == src].groupby("region_idx")["label"].agg(lambda s: ", ".join(sorted(set(s))))
        ).fillna("")
cols = (["contig", "region_idx", "start", "end", "region_len", "region_type", "hgt1_rank",
         "region_taxon"] + RANKS + [c for c in ["CARD", "RIP-DB"] if c in members.columns])
st.dataframe(members[cols].rename(columns={r: f"contig {RANK_NAMES[r].lower()}" for r in RANKS}),
             hide_index=True, width="stretch")


def _show_group():
    st.session_state["f_groups"] = [gid]


if st.button(f"Show group {gid}'s contigs in the browser", on_click=_show_group):
    st.switch_page("views/browser.py")
