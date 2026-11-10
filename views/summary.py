import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from hgtviz import RANK_NAMES, RANKS
from hgtviz.genemap import CARD, RIP, SEQ_BLUE

store = st.session_state["store"]
filtered: pd.DataFrame = st.session_state["filtered"]
spec = st.session_state["spec"]

st.title("Summary")
st.caption("Everything on this page describes the contigs that pass the sidebar filters.")
if filtered.empty:
    st.warning("No contigs match the current filters.")
    st.stop()

names = set(filtered["contig"])
R = store.subset("regions", filtered)
H1 = R[R["hgt1_rank"].notna()]
E = store.subset("link_ends", filtered, hgt2_only=True,
                 cols=["contig", "partner"] + [f"p_{r}" for r in RANKS])
A = store.subset("annotations", filtered) if store.has("annotations") else None

m = st.columns(5)
m[0].metric("Contigs", f"{len(filtered):,}")
m[1].metric("HGT1 regions", f"{len(H1):,}")
lo = np.minimum(E["_cid"].to_numpy(), E["_pid"].to_numpy()).astype(np.int64)
hi = np.maximum(E["_cid"].to_numpy(), E["_pid"].to_numpy()).astype(np.int64)
n_pairs = len(np.unique(lo * (len(store.ids()) + 1) + hi))
m[2].metric("HGT2 contig pairs", f"{n_pairs:,}",
            help="Distinct contig pairs with at least one contig in the filtered set.")
m[3].metric("Donor–recipient pairs",
            f"{store.dr[store.dr['recipient'].isin(names) | store.dr['donor'].isin(names)].shape[0]:,}")
m[4].metric("CARD hits in regions",
            "—" if A is None else f"{int(((A['source'] == 'card') & A['in_region']).sum()):,}")


def _layout(fig, h):
    fig.update_layout(height=h, margin=dict(l=10, r=10, t=30, b=10),
                      plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)")
    return fig


# ── HGT call ranks ───────────────────────────────────────────────────────────
st.subheader("Rank of the HGT calls")
rank_counts = pd.DataFrame({
    "HGT1 (contigs)": filtered["hgt1_rank"].value_counts(),
    "HGT2 (contigs)": filtered["hgt2_rank"].value_counts(),
}).reindex(RANKS).fillna(0).astype(int)
fig = go.Figure()
for col, color in zip(rank_counts.columns, [CARD, RIP]):
    fig.add_bar(x=[RANK_NAMES[r] for r in rank_counts.index], y=rank_counts[col], name=col,
                marker_color=color, hovertemplate="%{x}: %{y} contigs<extra>" + col + "</extra>")
fig.update_layout(barmode="group", bargap=0.3, bargroupgap=0.08, legend=dict(orientation="h", y=1.12),
                  yaxis_title="Contigs", xaxis_title="Coarsest rank at which the call was made")
st.plotly_chart(_layout(fig, 320), width="stretch")

# ── taxon matrices ───────────────────────────────────────────────────────────
rank = st.selectbox("Rank for the matrices", RANKS, index=RANKS.index(spec.tax_rank),
                    format_func=RANK_NAMES.get, key="sum_rank")
top_n = st.slider("Top taxa per axis", 5, 40, 15, key="sum_top")


def heatmap(df, x, y, xlab, ylab, symmetric=False):
    if df.empty:
        st.caption("Nothing to show at this rank.")
        return
    counts = df.groupby([y, x]).size().rename("n").reset_index()
    ty = counts.groupby(y)["n"].sum().nlargest(top_n).index
    tx = counts.groupby(x)["n"].sum().nlargest(top_n).index
    if symmetric:
        tx = ty = pd.Index(sorted(set(tx) | set(ty), key=lambda t: -counts.loc[
            (counts[x] == t) | (counts[y] == t), "n"].sum()))[:top_n]
    mat = (counts[counts[y].isin(ty) & counts[x].isin(tx)]
           .pivot(index=y, columns=x, values="n").reindex(index=ty, columns=tx).fillna(0))
    z = mat.where(mat > 0)  # empty cells stay blank instead of the lightest blue
    fig = go.Figure(go.Heatmap(z=z.values, x=list(mat.columns), y=list(mat.index),
                               colorscale=SEQ_BLUE, xgap=2, ygap=2,
                               hovertemplate=f"{ylab}: %{{y}}<br>{xlab}: %{{x}}<br>%{{z}}<extra></extra>"))
    fig.update_layout(xaxis_title=xlab, yaxis_title=ylab, yaxis_autorange="reversed",
                      xaxis_tickangle=-40)
    st.plotly_chart(_layout(fig, max(360, 26 * len(mat))), width="stretch")
    st.caption(f"{int(mat.values.sum()):,} of {len(df):,} shown (top {top_n} per axis).")


left, right = st.columns(2)
with left:
    st.markdown("**HGT1: region taxon → host contig taxon**")
    st.caption("One count per HGT1 region. Rows are the region's (putative donor) lineage.")
    h = pd.DataFrame({f"r_{rank}": H1[f"r_{rank}"].to_numpy(),
                      rank: store.contigs[rank].to_numpy()[H1["_cid"].to_numpy()]}).dropna()
    h = h[h[f"r_{rank}"] != h[rank]]
    heatmap(h, rank, f"r_{rank}", "Host contig", "Region")
with right:
    st.markdown("**HGT2: contig taxon ↔ partner taxon**")
    st.caption("One count per distinct contig pair; HGT2 has no direction, so the matrix is symmetric.")
    pairs = E[["_cid", "_pid", f"p_{rank}"]].drop_duplicates(["_cid", "_pid"])
    pairs = pd.DataFrame({f"p_{rank}": pairs[f"p_{rank}"].to_numpy(),
                          rank: store.contigs[rank].to_numpy()[pairs["_cid"].to_numpy()]}).dropna()
    pairs = pairs[pairs[rank] != pairs[f"p_{rank}"]]
    heatmap(pairs, f"p_{rank}", rank, "Partner", "Contig", symmetric=True)

# ── function ─────────────────────────────────────────────────────────────────
st.subheader("Function")
if A is None or A.empty:
    st.caption("No CARD / RIP-DB annotations on these contigs.")
    st.stop()


def where_bars(df, key, title, explode=False):
    if df.empty:
        st.caption(f"No {title.lower()}.")
        return
    d = df[["contig", key, "in_region"]].dropna()
    if explode:
        d = d.assign(**{key: d[key].str.split(";")}).explode(key)
    # one count per contig and label; "in a region" wins if a contig has both
    d = d.sort_values("in_region", ascending=False).drop_duplicates(["contig", key])
    t = (d.groupby([key, "in_region"]).size().unstack(fill_value=0)
         .reindex(columns=[True, False], fill_value=0))
    t = t.loc[t.sum(axis=1).sort_values().index].tail(25)
    fig = go.Figure()
    for flag, name, color in [(True, "inside a kairos region", CARD), (False, "elsewhere only", RIP)]:
        fig.add_bar(y=t.index, x=t[flag], name=name, orientation="h", marker_color=color,
                    hovertemplate="%{y}: %{x} contigs<extra>" + name + "</extra>")
    fig.update_layout(barmode="stack", bargap=0.25, xaxis_title="Contigs",
                      legend=dict(orientation="h", yanchor="top", y=-0.18, x=0))
    st.markdown(f"**{title}**")
    st.plotly_chart(_layout(fig, max(300, 22 * len(t) + 120)), width="stretch")


card = A[A["source"] == "card"]
c1, c2, c3 = st.columns(3)
with c1:
    where_bars(card, "label", "CARD genes")
with c2:
    if card["category"].notna().any():
        where_bars(card, "category", "CARD drug classes", explode=True)
    else:
        st.caption("No drug classes (build the store with --card-index).")
with c3:
    where_bars(A[A["source"] == "ripdb"], "label", "RIP-DB families")
