"""Taxonomy helpers shared by the build step and the app."""
import re

import numpy as np
import pandas as pd

from . import RANKS

MMSEQS_COLS = ["id", "taxid", "tax_level", "tax_name", "fragments",
               "returned", "agreed", "percent_agreed", "lineage"]
_PREFIX = re.compile(r"^([a-z])_(.*)$")


def parse_lineage(lineage) -> dict:
    """'d_Bacteria;p_Bacillota;...' -> {'d': 'Bacteria', 'p': 'Bacillota', ...}.

    Only the leading rank prefix of each element is removed, so GTDB names that
    contain '<letter>_' (e.g. 'Pseudomonas_A') are kept intact. Elements with a
    prefix that is not one of the seven ranks (e.g. '-_root') are skipped.
    """
    out = {}
    if not isinstance(lineage, str):
        return out
    for part in lineage.split(";"):
        m = _PREFIX.match(part.strip())
        if m and m.group(1) in RANKS and m.group(2):
            out[m.group(1)] = m.group(2)
    return out


def read_mmseqs_taxonomy(path, ids=None, id_transform=None) -> pd.DataFrame:
    """Read an mmseqs `createtsv` taxonomy table (no header, --tax-lineage 1).

    ids: optional collection; only rows whose (transformed) id is in it are parsed.
    Unclassified rows (taxid <= 1) get empty ranks.
    """
    df = pd.read_csv(path, sep="\t", header=None, names=MMSEQS_COLS,
                     dtype={"id": str, "lineage": str}, na_filter=False)
    if id_transform is not None:
        df["id"] = id_transform(df["id"])
    if ids is not None:
        df = df[df["id"].isin(set(ids))]
    df = df.drop_duplicates("id").reset_index(drop=True)
    df["taxid"] = pd.to_numeric(df["taxid"], errors="coerce").fillna(0).astype(int)
    for col in ["fragments", "returned", "agreed", "percent_agreed"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    parsed = pd.DataFrame([parse_lineage(x) for x in df["lineage"]],
                          columns=RANKS, index=df.index)
    parsed[df["taxid"] <= 1] = np.nan
    return pd.concat([df.drop(columns="lineage"), parsed], axis=1)


def discordant_rank(a: pd.DataFrame, b: pd.DataFrame, a_fmt="{}", b_fmt="{}") -> pd.Series:
    """Highest rank at which two lineages are both named and differ, given that
    every higher rank is named and equal (kairos `hgt_levels` logic).
    NaN when the lineages never become discordant."""
    out = pd.Series(np.nan, index=a.index, dtype=object)
    still_equal = pd.Series(True, index=a.index)
    for r in RANKS:
        x, y = a[a_fmt.format(r)], b[b_fmt.format(r)]
        both = x.notna() & y.notna()
        hit = still_equal & both & (x != y) & out.isna()
        out[hit] = r
        still_equal &= both & (x == y)
    return out


def rank_index(rank) -> float:
    return RANKS.index(rank) if rank in RANKS else np.inf


def coarsest(ranks) -> str | float:
    ranks = [r for r in ranks if r in RANKS]
    return min(ranks, key=RANKS.index) if ranks else np.nan
