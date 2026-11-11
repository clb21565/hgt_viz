"""Side-by-side lineages of the taxa involved in one contig's HGT calls: the host
contig, each HGT1 region (its own classification) and each distinct HGT2 partner
lineage, aligned by rank so the point of divergence is easy to see."""
import pandas as pd

from . import HGT1_LABEL, HGT2_LABEL, RANK_NAMES, RANKS

ORDER = {r: i for i, r in enumerate(RANKS)}
# cell styles (work on light and dark backgrounds)
SAME = "color: #8a8984"                                              # agrees with the host
SPLIT = "background-color: rgba(235,104,52,0.38); font-weight: 600"  # first rank that differs
BELOW = "background-color: rgba(235,104,52,0.14)"                    # differs below the split


def _rank_name(r):
    return RANK_NAMES.get(r, "") if isinstance(r, str) else ""


def hgt_lineages(host: pd.Series, regions: pd.DataFrame, links: pd.DataFrame | None,
                 use=("hgt1", "hgt2"), max_partner_rows: int = 30,
                 full: bool = False) -> tuple[pd.DataFrame, int]:
    """One row for the host, one per HGT1 region, one per distinct HGT2 partner
    lineage (partners grouped, most common first within each rank). Unless `full`,
    partner lineages are cut below the rank where they differ from the host, so
    partners are grouped by the diverging taxon. Columns:
    evidence, region(s), contigs, differs at, then one column per rank name.
    Returns the table and the number of partner lineages left out."""
    names = [RANK_NAMES[r] for r in RANKS]

    def ranks_of(row, prefix):
        return {RANK_NAMES[r]: row.get(f"{prefix}{r}") if pd.notna(row.get(f"{prefix}{r}")) else ""
                for r in RANKS}

    rows = [{"evidence": "this contig (host)", "region(s)": "", "contigs": None,
             "differs at": "", **ranks_of(host, "")}]
    if "hgt1" in use:
        h1 = regions[regions["hgt1_rank"].notna()].sort_values("start")
        for _, r in h1.iterrows():
            rows.append({"evidence": HGT1_LABEL, "region(s)": r["region_idx"], "contigs": None,
                         "differs at": _rank_name(r["hgt1_rank"]), **ranks_of(r, "r_")})
    hidden = 0
    if "hgt2" in use and links is not None and len(links):
        L = links[links["is_hgt2"]]
        if len(L):
            pcols = [f"p_{r}" for r in RANKS]
            L = L.assign(**{c: L[c].astype("object").where(L[c].notna(), "") for c in pcols})
            if not full:   # blank every rank below the split
                cut = L["discordant_rank"].map(lambda x: ORDER.get(x, 99))
                for i, c in enumerate(pcols):
                    L[c] = L[c].where(cut >= i, "")
            g = (L.groupby(pcols, sort=False)
                 .agg(contigs=("partner", "nunique"),
                      regions=("region", lambda s: sorted(set(s))),
                      rank=("discordant_rank", lambda s: min(s, key=lambda x: ORDER.get(x, 99))))
                 .reset_index())
            g["_o"] = g["rank"].map(lambda x: ORDER.get(x, 99))
            g = g.sort_values(["_o", "contigs"], ascending=[True, False])
            hidden = max(0, len(g) - max_partner_rows)
            for _, r in g.head(max_partner_rows).iterrows():
                regs = r["regions"]
                rows.append({"evidence": HGT2_LABEL,
                             "region(s)": ", ".join(regs[:3]) + (f" +{len(regs) - 3}" if len(regs) > 3 else ""),
                             "contigs": int(r["contigs"]), "differs at": _rank_name(r["rank"]),
                             **{RANK_NAMES[x]: r[f"p_{x}"] for x in RANKS}})
    df = pd.DataFrame(rows, columns=["evidence", "region(s)", "contigs", "differs at"] + names)
    df["contigs"] = df["contigs"].map(lambda x: "" if x is None or pd.isna(x) else str(int(x)))
    return df, hidden


def style(df: pd.DataFrame):
    """Grey where a lineage agrees with the host, strong highlight at the rank where
    it first differs, light highlight below that."""
    names = [RANK_NAMES[r] for r in RANKS]
    host = df.iloc[0]

    def row_style(row):
        out = {c: "" for c in df.columns}
        if row.name == df.index[0]:
            return pd.Series(out)
        split = row["differs at"]
        seen_split = False
        for c in names:
            if c == split:
                seen_split = True
                out[c] = SPLIT
            elif not seen_split:
                out[c] = SAME if row[c] == host[c] else ""
            else:
                out[c] = BELOW if row[c] and row[c] != host[c] else SAME
        return pd.Series(out)

    return df.style.apply(row_style, axis=1).format(na_rep="")
