"""Contig filters. Every filter narrows the set of base contigs; empty settings
are ignored. Kept free of Streamlit so it can be unit-tested."""
from dataclasses import astuple, dataclass, field

import re

import pandas as pd

from . import RANKS
from .store import Store

FUNCTION_FIELDS = ["card_genes", "drug_classes", "rip_families", "cog_letters", "eggnog_text",
                   "function_scope", "function_hgt_ranks"]
SCOPES = {"contig": "anywhere on the contig", "region": "inside a kairos region",
          "hgt": "inside a region with an HGT call"}


@dataclass
class FilterSpec:
    search: str = ""                 # substring of the contig ID
    require: list = field(default_factory=list)  # any of "HGT1", "HGT2", "Donor–recipient", "In situ"
    hgt1_ranks: list = field(default_factory=list)
    hgt2_ranks: list = field(default_factory=list)
    tax_rank: str = "p"
    contig_taxa: list = field(default_factory=list)    # host contig taxon at tax_rank
    region_taxa: list = field(default_factory=list)    # HGT1 region taxon at tax_rank
    partner_taxa: list = field(default_factory=list)   # HGT2 partner taxon at tax_rank
    host_lineage: str = ""           # text in the host contig's full lineage
    other_lineage: str = ""          # text in an HGT1 region's or HGT2 partner's full lineage
    other_lineage_in: str = "either"  # "either", "hgt1" or "hgt2"
    groups: list = field(default_factory=list)
    min_group_size: int = 0
    card_genes: list = field(default_factory=list)
    drug_classes: list = field(default_factory=list)
    rip_families: list = field(default_factory=list)
    cog_letters: list = field(default_factory=list)   # COG categories (any of)
    eggnog_text: str = ""            # substring of eggNOG name/description/KO/Pfam/CAZy/EC
    function_scope: str = "contig"   # "contig", "region" (inside a kairos region) or "hgt"
                                     # (inside a region with an HGT1 call or HGT2 link)
    function_hgt_ranks: list = field(default_factory=list)  # scope "hgt": coarsest HGT rank in
    min_length: int = 0
    min_partners: int = 0

    def active(self) -> list[str]:
        """Human-readable list of the filters that are set."""
        out = []
        for name in ["search", "require", "hgt1_ranks", "hgt2_ranks", "contig_taxa", "region_taxa",
                     "partner_taxa", "host_lineage", "other_lineage", "groups", "card_genes", "drug_classes", "rip_families",
                     "cog_letters", "eggnog_text"]:
            if getattr(self, name):
                out.append(name)
        for name in ["min_group_size", "min_length", "min_partners"]:
            if getattr(self, name):
                out.append(name)
        return out


def _contigs_with(df: pd.DataFrame, mask) -> set:
    return set(df.loc[mask, "contig"])


def hgt_regions(store: Store, spec: FilterSpec) -> frozenset | None:
    """For scope 'hgt': regions with an HGT1 call or HGT2 link (coarsest rank among
    function_hgt_ranks, if set). None for the other scopes."""
    if spec.function_scope != "hgt":
        return None
    def run():
        h = store.region_hgt()
        m = h["hgt_rank"].notna()
        if spec.function_hgt_ranks:
            m &= h["hgt_rank"].isin(spec.function_hgt_ranks)
        return frozenset(h.index[m])
    return store.memo(("hgt_regions", tuple(spec.function_hgt_ranks)), run, maxsize=8)


def annotation_hits(store: Store, spec: FilterSpec) -> pd.DataFrame:
    a = store.annotations
    if a is None:
        return pd.DataFrame(columns=["contig", "source", "label", "category", "in_region",
                                     "region_idx", "family", "orf_id"])
    if spec.function_scope == "region":
        return a[a["in_region"]]
    if spec.function_scope == "hgt":
        return a[a["region_idx"].isin(hgt_regions(store, spec))]
    return a


def lineage_terms(text: str) -> str | None:
    """Comma-separated terms (any of them) as one case-insensitive regex."""
    terms = [t.strip() for t in text.split(",") if t.strip()]
    return "|".join(re.escape(t) for t in terms) or None


def lineage(store: Store, name: str, prefix: str = "") -> pd.Series:
    """Full lineage of every row of a table ('d; p; c; o; f; g; s', blanks dropped),
    built once per table and column prefix."""
    def run():
        df = store.link_ends if name == "link_ends" else getattr(store, name)
        cols = [f"{prefix}{r}" for r in RANKS if f"{prefix}{r}" in df.columns]
        x = df[cols].astype("object").where(df[cols].notna(), "").astype(str)
        out = x[cols[0]] if cols else pd.Series("", index=df.index)
        for c in cols[1:]:
            out = out + "; " + x[c]
        return out.str.replace(r"(; )+$", "", regex=True) + ";"   # 'Bacillota;' = exact name
    return store.memo(("lineage", name, prefix), run, maxsize=8, bucket="lineage")


def lineage_matches(store: Store, name: str, prefix: str, text: str) -> pd.Series:
    pat = lineage_terms(text)
    lin = lineage(store, name, prefix)
    return lin.str.contains(pat, case=False, regex=True) if pat else pd.Series(True, index=lin.index)


def _key(spec: FilterSpec):
    return tuple(tuple(v) if isinstance(v, list) else v for v in astuple(spec))


def apply(store: Store, spec: FilterSpec) -> pd.DataFrame:
    """Return the rows of store.contigs that pass every filter (memoised per spec)."""
    return store.memo(("apply", _key(spec)), lambda: _apply(store, spec), maxsize=16, bucket="apply")


def _apply(store: Store, spec: FilterSpec) -> pd.DataFrame:
    c = store.contigs
    keep = pd.Series(True, index=c.index)

    if spec.search:
        keep &= c["contig"].str.contains(spec.search, case=False, regex=False)
    req = set(spec.require)
    if "HGT1" in req:
        keep &= c["hgt1"]
    if "HGT2" in req:
        keep &= c["hgt2"]
    if "Donor–recipient" in req:
        keep &= c["dr_role"] != ""
    if "In situ" in req and "insitu_as_recipient" in c.columns:
        keep &= (c["insitu_as_recipient"] + c["insitu_as_donor"]) > 0
    if spec.min_length:
        keep &= c["length"] >= spec.min_length
    if spec.min_partners:
        keep &= c["n_hgt2_partners"] >= spec.min_partners

    r = spec.tax_rank
    if spec.contig_taxa:
        keep &= c[r].isin(spec.contig_taxa)
    if lineage_terms(spec.host_lineage):
        keep &= lineage_matches(store, "contigs", "", spec.host_lineage)

    sets = []
    R = store.regions
    if spec.hgt1_ranks:
        sets.append(_contigs_with(R, R["hgt1_rank"].isin(spec.hgt1_ranks)))
    if spec.region_taxa:
        sets.append(_contigs_with(R, R["hgt1_rank"].notna() & R[f"r_{r}"].isin(spec.region_taxa)))
    if spec.groups:
        sets.append(_contigs_with(R, R["group_id"].isin(spec.groups)))
    if spec.min_group_size:
        sets.append(_contigs_with(R, R["group_size"] >= spec.min_group_size))

    if spec.hgt2_ranks or spec.partner_taxa:
        E = store.link_ends
        E = E[E["is_hgt2"]]
        if spec.hgt2_ranks:
            sets.append(_contigs_with(E, E["discordant_rank"].isin(spec.hgt2_ranks)))
        if spec.partner_taxa:
            sets.append(_contigs_with(E, E[f"p_{r}"].isin(spec.partner_taxa)))

    if lineage_terms(spec.other_lineage):
        hits = set()
        if spec.other_lineage_in in ("either", "hgt1"):
            m = R["hgt1_rank"].notna() & lineage_matches(store, "regions", "r_", spec.other_lineage)
            hits |= _contigs_with(R, m)
        if spec.other_lineage_in in ("either", "hgt2"):
            E = store.link_ends
            m = E["is_hgt2"] & lineage_matches(store, "link_ends", "p_", spec.other_lineage)
            hits |= _contigs_with(E, m)
        sets.append(hits)

    if spec.card_genes or spec.drug_classes or spec.rip_families:
        A = annotation_hits(store, spec)
        if spec.card_genes:
            sets.append(_contigs_with(A, (A["source"] == "card") & A["label"].isin(spec.card_genes)))
        if spec.drug_classes:
            dc = A[["contig", "category"]].assign(category=A["category"].fillna("").str.split(";"))
            dc = dc[(A["source"] == "card").to_numpy()].explode("category")
            m = dc["category"].isin(spec.drug_classes)
            sets.append(_contigs_with(dc, m))
        if spec.rip_families:
            sets.append(_contigs_with(A, (A["source"] == "ripdb") & A["label"].isin(spec.rip_families)))

    if spec.cog_letters and spec.function_scope != "hgt" and store.cog_counts is not None:
        cc = store.cog_counts
        m = cc["cog"].isin(spec.cog_letters)
        if spec.function_scope == "region":
            m &= cc["in_region"]
        sets.append(_contigs_with(cc, m))
    elif spec.cog_letters and store.has("eggnog"):    # per gene: needs the region of each gene
        sets.append(store.eggnog_search("", tuple(spec.cog_letters), regions=hgt_regions(store, spec)))
    if spec.eggnog_text.strip():
        sets.append(store.eggnog_search(spec.eggnog_text.strip(), (), spec.function_scope == "region",
                                        hgt_regions(store, spec)))

    for s in sets:
        keep &= c["contig"].isin(s)
    return c[keep]


def options(store: Store, rank: str) -> dict:
    """Choices for the multiselects (memoised per rank)."""
    return store.memo(("options", rank), lambda: _options(store, rank))


def _options(store: Store, rank: str) -> dict:
    c, R = store.contigs, store.regions
    E = store.link_ends
    E = E[E["is_hgt2"]]
    a = store.annotations
    out = {
        "contig_taxa": sorted(c[rank].dropna().unique()),
        "region_taxa": sorted(R.loc[R["hgt1_rank"].notna(), f"r_{rank}"].dropna().unique()),
        "partner_taxa": sorted(E[f"p_{rank}"].dropna().unique()),
        "groups": sorted(int(g) for g in R["group_id"].dropna().unique()),
        "card_genes": [], "drug_classes": [], "rip_families": [], "cog_letters": [],
    }
    if store.cog_counts is not None:
        from .cog import COG
        present = set(store.cog_counts["cog"].unique())
        out["cog_letters"] = [x for x in COG if x in present]
    if a is not None:
        card = a[a["source"] == "card"]
        out["card_genes"] = sorted(card["label"].dropna().unique())
        out["drug_classes"] = sorted({x for s in card["category"].dropna() for x in s.split(";") if x})
        out["rip_families"] = sorted(a.loc[a["source"] == "ripdb", "label"].dropna().unique())
    return out


def function_active(spec: FilterSpec) -> bool:
    return bool(spec.card_genes or spec.drug_classes or spec.rip_families or spec.cog_letters
                or spec.eggnog_text.strip())


def function_hits(store: Store, spec: FilterSpec, contigs: set | None = None,
                  limit: int = 5000) -> tuple[pd.DataFrame, dict]:
    """Genes that match any of the function criteria in `spec` (on `contigs`, if
    given). One row per gene and matching source: orf_id, contig, region_idx,
    in_region, source, name, detail. Also returns the total number of matching
    genes per source (the table is capped at `limit` rows per source)."""
    frames, totals = [], {}
    region = spec.function_scope == "region"
    A = annotation_hits(store, spec)
    if contigs is not None:
        A = A[A["contig"].isin(contigs)]
    for src, labels, name in [("card", spec.card_genes, "CARD gene"),
                              ("ripdb", spec.rip_families, "RIP-DB family")]:
        if labels:
            a = A[(A["source"] == src) & A["label"].isin(labels)]
            totals[name] = len(a)
            frames.append(a.head(limit).assign(source=name, name=a["label"],
                                               detail=a["family"].fillna("")))
    if spec.drug_classes:
        a = A[A["source"] == "card"]
        m = a["category"].fillna("").str.split(";").map(lambda xs: bool(set(xs) & set(spec.drug_classes)))
        a = a[m.astype(bool)]
        totals["CARD drug class"] = len(a)
        frames.append(a.head(limit).assign(source="CARD drug class", name=a["label"],
                                           detail=a["category"].fillna("")))
    text = spec.eggnog_text.strip()
    if store.has("eggnog") and (text or spec.cog_letters):
        e, n = store.eggnog_rows(text, tuple(spec.cog_letters), region, contigs, limit,
                                 hgt_regions(store, spec))
        totals["eggNOG"] = n
        if len(e):
            desc = e["description"].fillna("")
            frames.append(e.assign(
                source="eggNOG", name=e["preferred_name"].fillna(""),
                detail=desc + e["cog_category"].map(lambda x: f" [COG {x}]" if isinstance(x, str) else "")))
    cols = ["orf_id", "contig", "region_idx", "in_region", "source", "name", "detail"]
    if not frames:
        out = pd.DataFrame(columns=cols)
    else:
        out = pd.concat([f[cols] for f in frames], ignore_index=True)
    h = store.region_hgt()
    for col in ["hgt1_rank", "hgt2_rank", "hgt2_partners", "hgt_rank"]:
        out[col] = out["region_idx"].map(h[col])
    out["hgt2_partners"] = out["hgt2_partners"].fillna(0).astype(int)
    return out, totals
