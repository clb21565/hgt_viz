"""Load a store written by hgtviz.build."""
import json
import os
import re
from collections import OrderedDict
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from . import RANKS

TABLES = ["contigs", "regions", "links", "dr", "groups", "genes", "annotations", "insitu",
          "cog_counts"]
LAZY = ["eggnog"]   # large; read per contig from Parquet, never loaded whole
EGGNOG_SEARCH_COLS = ["preferred_name", "description", "kegg_ko", "pfams", "cazy", "ec"]


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
    cog_counts: pd.DataFrame | None = None
    _ends: pd.DataFrame | None = field(default=None, repr=False)
    _index: dict = field(default_factory=dict, repr=False)
    caches: dict = field(default_factory=dict, repr=False)

    def memo(self, key, fn, maxsize=64, bucket="small"):
        """LRU shared by all sessions (the store is read-only). Large results go in
        their own bucket with a small maxsize to bound memory."""
        cache = self.caches.setdefault(bucket, OrderedDict())
        if key in cache:
            cache.move_to_end(key)
            return cache[key]
        val = fn()
        cache[key] = val
        if len(cache) > maxsize:
            cache.popitem(last=False)
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

    def ids(self) -> pd.Index:
        """Every contig name in the store. Base contigs come first, in store.contigs
        order, so a base contig's code equals its row position in store.contigs."""
        if "ids" not in self._index:
            base = pd.Index(self.contigs["contig"])
            other = pd.Index(pd.unique(pd.concat(
                [self.links["contig_a"], self.links["contig_b"], self.regions["contig"]],
                ignore_index=True)))
            self._index["ids"] = base.append(other.difference(base))
        return self._index["ids"]

    def codes(self, name, col="contig") -> np.ndarray:
        """Integer contig code of each row of table `name` (cached)."""
        key = ("codes", name, col)
        if key not in self._index:
            df = self.link_ends if name == "link_ends" else getattr(self, name)
            self._index[key] = self.ids().get_indexer(df[col])
        return self._index[key]

    def subset(self, name, filtered: pd.DataFrame, hgt2_only=False, cols=None) -> pd.DataFrame:
        """Rows of `name` (regions, link_ends, annotations) on the filtered contigs, with
        integer codes `_cid` (and `_pid` for the partner of a link end). A boolean
        lookup on codes instead of a string isin. Only the row positions are memoised
        (per filtered set); pass `cols` to copy just the columns you need.
        `filtered` must be a slice of store.contigs (its index = row positions)."""
        sig = hash(filtered.index.to_numpy().tobytes())
        pos = self.memo(("subset", name, hgt2_only, len(filtered), sig),
                        lambda: self._positions(name, filtered, hgt2_only), maxsize=32,
                        bucket="subsets")
        df = self.link_ends if name == "link_ends" else getattr(self, name)
        out = (df if cols is None else df[[c for c in cols if c in df.columns]]).iloc[pos].copy()
        out["_cid"] = self.codes(name)[pos]
        if name == "link_ends":
            out["_pid"] = self.codes(name, "partner")[pos]
        return out

    def _positions(self, name, filtered, hgt2_only) -> np.ndarray:
        df = self.link_ends if name == "link_ends" else getattr(self, name)
        cid = self.codes(name)
        if len(filtered) == len(self.contigs):          # no filter: every base contig
            rows = cid >= 0
            rows &= cid < len(self.contigs)
        else:
            mask = np.zeros(len(self.ids()) + 1, dtype=bool)   # last slot: code -1
            mask[filtered.index.to_numpy()] = True
            rows = mask[cid]
        if hgt2_only:
            rows &= df["is_hgt2"].to_numpy(dtype=bool)
        return np.flatnonzero(rows)

    # ── eggNOG (lazy) ──────────────────────────────────────────────────────────
    @property
    def eggnog_path(self) -> str | None:
        p = os.path.join(self.path, "eggnog.parquet")
        return p if os.path.exists(p) else None

    def eggnog_for(self, contig: str) -> pd.DataFrame | None:
        """eggNOG rows of one contig. The file is sorted by contig with small row
        groups, so the min/max statistics let pyarrow read one or two of them."""
        if self.eggnog_path is None:
            return None
        return self.memo(("eggnog", contig), lambda: pq.read_table(
            self.eggnog_path, filters=[("contig", "==", contig)]).to_pandas(),
            maxsize=256, bucket="eggnog")

    @staticmethod
    def _eggnog_mask(t, text, cog_letters, region_only, regions):
        m = pa.array([True] * t.num_rows)
        if text:
            tm = None
            for c in EGGNOG_SEARCH_COLS:
                x = pc.fill_null(pc.match_substring(t[c], text, ignore_case=True), False)
                tm = x if tm is None else pc.or_(tm, x)
            m = pc.and_(m, tm)
        if cog_letters:
            cm = None
            for L in cog_letters:
                x = pc.fill_null(pc.match_substring(t["cog_category"], L), False)
                cm = x if cm is None else pc.or_(cm, x)
            m = pc.and_(m, cm)
        if region_only:
            m = pc.and_(m, pc.fill_null(t["in_region"], False))
        if regions is not None:
            vs = pa.array(sorted(regions), pa.string())
            m = pc.and_(m, pc.fill_null(pc.is_in(t["region_idx"], value_set=vs), False))
        return m

    def eggnog_search(self, text: str = "", cog_letters=(), region_only=False,
                      regions: frozenset | None = None) -> set:
        """Contigs with a gene whose eggNOG name, description, KO, Pfam, CAZy or EC
        contains `text` (case-insensitive) and/or has any of the COG letters,
        optionally only inside kairos regions or inside the given regions."""
        def run():
            hits = set()
            cols = ["contig", "in_region", "region_idx", "cog_category"] + EGGNOG_SEARCH_COLS
            pf = pq.ParquetFile(self.eggnog_path)
            for i in range(pf.num_row_groups):
                t = pf.read_row_group(i, columns=list(dict.fromkeys(cols)))
                m = self._eggnog_mask(t, text, cog_letters, region_only, regions)
                hits.update(pc.unique(pc.filter(t["contig"], m)).to_pylist())
            return hits
        if self.eggnog_path is None or not (text or cog_letters):
            return set()
        key = ("eggnog_search", text.lower(), tuple(cog_letters), region_only,
               None if regions is None else hash(regions))
        return self.memo(key, run, maxsize=32, bucket="eggnog_search")

    def eggnog_rows(self, text: str = "", cog_letters=(), region_only=False,
                    contigs: set | None = None, limit: int = 5000,
                    regions: frozenset | None = None) -> tuple[pd.DataFrame, int]:
        """Genes whose eggNOG annotation contains `text` and/or has any of the COG
        letters, optionally only on `contigs`, inside kairos regions or inside the
        given regions. Returns at most `limit` rows and the total number of matches."""
        def run():
            cols = ["orf_id", "contig", "region_idx", "in_region", "cog_category"] + EGGNOG_SEARCH_COLS
            pf = pq.ParquetFile(self.eggnog_path)
            vs = pa.array(sorted(contigs), pa.string()) if contigs is not None else None
            parts, total = [], 0
            for i in range(pf.num_row_groups):
                t = pf.read_row_group(i, columns=list(dict.fromkeys(cols)))
                m = self._eggnog_mask(t, text, cog_letters, region_only, regions)
                if vs is not None:
                    m = pc.and_(m, pc.is_in(t["contig"], value_set=vs))
                n = pc.sum(m).as_py() or 0
                total += n
                if n and sum(len(x) for x in parts) < limit:
                    parts.append(pc.filter(t, m).to_pandas())
            out = pd.concat(parts, ignore_index=True).head(limit) if parts else pd.DataFrame(columns=cols)
            return out, total
        if self.eggnog_path is None or not (text or cog_letters):
            return pd.DataFrame(), 0
        key = ("eggnog_rows", text.lower(), tuple(cog_letters), region_only,
               None if contigs is None else hash(frozenset(contigs)), limit,
               None if regions is None else hash(regions))
        return self.memo(key, run, maxsize=8, bucket="eggnog_rows")

    def region_hgt(self) -> pd.DataFrame:
        """HGT evidence per kairos region (index region_idx): hgt1_rank, hgt2_rank
        (coarsest rank over its HGT2 links), hgt2_partners and hgt_rank (coarsest
        of the two; missing = region without an HGT call)."""
        def run():
            R = self.regions.set_index("region_idx")
            E = self.link_ends
            E = E[E["is_hgt2"]]
            order = {r: i for i, r in enumerate(RANKS)}
            ri = E["discordant_rank"].map(order)
            h2 = (E.assign(_ri=ri).groupby("region", observed=True)
                  .agg(_ri=("_ri", "min"), hgt2_partners=("partner", "nunique")))
            out = pd.DataFrame(index=R.index)
            out["hgt1_rank"] = R["hgt1_rank"]
            out["hgt2_rank"] = h2["_ri"].reindex(out.index).map(
                lambda i: RANKS[int(i)] if pd.notna(i) else None)
            out["hgt2_partners"] = h2["hgt2_partners"].reindex(out.index).fillna(0).astype(int)
            both = pd.concat([out["hgt1_rank"].map(order), out["hgt2_rank"].map(order)], axis=1).min(axis=1)
            out["hgt_rank"] = both.map(lambda i: RANKS[int(i)] if pd.notna(i) else None)
            return out
        return self.memo(("region_hgt",), run, maxsize=1, bucket="region_hgt")

    def has(self, name) -> bool:
        if name == "eggnog":
            return self.eggnog_path is not None
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
