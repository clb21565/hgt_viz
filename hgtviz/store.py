"""Load a store written by hgtviz.build."""
import json
import os
import re
from collections import OrderedDict
from dataclasses import dataclass, field

import pandas as pd

from . import RANKS

TABLES = ["contigs", "regions", "links", "dr", "groups", "genes", "annotations", "insitu"]


class StoreError(Exception):
    pass


@dataclass
class Store:
    path: str
    meta: dict
    contigs: pd.DataFrame
    regions: pd.DataFrame
    links: pd.DataFrame
    dr: pd.DataFrame
    groups: pd.DataFrame | None = None
    genes: pd.DataFrame | None = None
    annotations: pd.DataFrame | None = None
    insitu: pd.DataFrame | None = None
    _ends: pd.DataFrame | None = field(default=None, repr=False)
    _index: dict = field(default_factory=dict, repr=False)
    cache: OrderedDict = field(default_factory=OrderedDict, repr=False)

    def memo(self, key, fn, maxsize=64):
        """Small LRU shared by all sessions (the store is read-only)."""
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        val = fn()
        self.cache[key] = val
        if len(self.cache) > maxsize:
            self.cache.popitem(last=False)
        return val

    def rows_for(self, name, contig, col="contig") -> pd.DataFrame | None:
        """Rows of table `name` whose `col` equals `contig`, via a positional index
        built once per table instead of a full scan per lookup."""
        df = self.link_ends if name == "link_ends" else getattr(self, name, None)
        if df is None:
            return None
        key = (name, col)
        if key not in self._index:
            self._index[key] = df.groupby(col, sort=False, observed=True).indices
        pos = self._index[key].get(contig)
        return df.iloc[pos] if pos is not None else df.iloc[:0]

    def has(self, name) -> bool:
        df = getattr(self, name, None)
        return df is not None and len(df) > 0

    @property
    def link_ends(self) -> pd.DataFrame:
        """Links in both orientations: one row per (contig, partner) end, with the
        partner's taxonomy as p_<rank> and this end's coordinates as start/end."""
        if self._ends is None:
            L = self.links
            parts = []
            for me, other in [("a", "b"), ("b", "a")]:
                d = pd.DataFrame({
                    "link": L.index,
                    "contig": L[f"contig_{me}"], "region": L[f"region_{me}"],
                    "start": L.get(f"start_{me}"), "end": L.get(f"end_{me}"),
                    "partner": L[f"contig_{other}"], "partner_region": L[f"region_{other}"],
                    "partner_start": L.get(f"start_{other}"), "partner_end": L.get(f"end_{other}"),
                })
                for r in RANKS:
                    d[f"p_{r}"] = L[f"{other}_{r}"]
                for c in ["strand", "identity", "shared_len", "n_shared_orfs",
                          "discordant_rank", "is_hgt2", "event_id"]:
                    if c in L.columns:
                        d[c] = L[c]
                parts.append(d)
            self._ends = pd.concat(parts, ignore_index=True)
        return self._ends


def is_store(path) -> bool:
    return os.path.exists(os.path.join(path, "meta.json"))


def load(path) -> Store:
    if not os.path.isdir(path):
        raise StoreError(f"Not a directory: {path}")
    meta_path = os.path.join(path, "meta.json")
    if not os.path.exists(meta_path):
        raise StoreError(f"No meta.json in {path}: build a store with `python -m hgtviz.build`")
    with open(meta_path) as fh:
        meta = json.load(fh)
    tables = {}
    for name in TABLES:
        p = os.path.join(path, f"{name}.parquet")
        tables[name] = pd.read_parquet(p) if os.path.exists(p) else None
    for req in ["contigs", "regions", "links", "dr"]:
        if tables[req] is None:
            raise StoreError(f"{path} is missing {req}.parquet; rebuild the store")
    tables["contigs"] = natural_sort(tables["contigs"])
    return Store(path=path, meta=meta, **tables)


def _natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def natural_sort(df: pd.DataFrame, col="contig") -> pd.DataFrame:
    """contig_2 before contig_10. Done once at load so filtered views keep the order."""
    order = sorted(range(len(df)), key=lambda i: _natural_key(str(df[col].iat[i])))
    return df.iloc[order].reset_index(drop=True)
