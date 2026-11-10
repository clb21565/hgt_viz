"""Contig filters. Every filter narrows the set of base contigs; empty settings
are ignored. Kept free of Streamlit so it can be unit-tested."""
from dataclasses import dataclass, field

import pandas as pd

from .store import Store


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
    groups: list = field(default_factory=list)
    min_group_size: int = 0
    card_genes: list = field(default_factory=list)
    drug_classes: list = field(default_factory=list)
    rip_families: list = field(default_factory=list)
    function_scope: str = "contig"   # "contig" or "region" (inside a kairos region)
    min_length: int = 0
    min_partners: int = 0

    def active(self) -> list[str]:
        """Human-readable list of the filters that are set."""
        out = []
        for name in ["search", "require", "hgt1_ranks", "hgt2_ranks", "contig_taxa", "region_taxa",
                     "partner_taxa", "groups", "card_genes", "drug_classes", "rip_families"]:
            if getattr(self, name):
                out.append(name)
        for name in ["min_group_size", "min_length", "min_partners"]:
            if getattr(self, name):
                out.append(name)
        return out


def _contigs_with(df: pd.DataFrame, mask) -> set:
    return set(df.loc[mask, "contig"])


def annotation_hits(store: Store, scope: str) -> pd.DataFrame:
    a = store.annotations
    if a is None:
        return pd.DataFrame(columns=["contig", "source", "label", "category", "in_region"])
    return a[a["in_region"]] if scope == "region" else a


def apply(store: Store, spec: FilterSpec) -> pd.DataFrame:
    """Return the rows of store.contigs that pass every filter."""
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

    if spec.card_genes or spec.drug_classes or spec.rip_families:
        A = annotation_hits(store, spec.function_scope)
        if spec.card_genes:
            sets.append(_contigs_with(A, (A["source"] == "card") & A["label"].isin(spec.card_genes)))
        if spec.drug_classes:
            dc = A["category"].fillna("").str.split(";")
            m = (A["source"] == "card") & dc.apply(lambda x: bool(set(x) & set(spec.drug_classes)))
            sets.append(_contigs_with(A, m))
        if spec.rip_families:
            sets.append(_contigs_with(A, (A["source"] == "ripdb") & A["label"].isin(spec.rip_families)))

    for s in sets:
        keep &= c["contig"].isin(s)
    return c[keep]


def options(store: Store, rank: str) -> dict:
    """Choices for the multiselects."""
    c, R = store.contigs, store.regions
    E = store.link_ends
    E = E[E["is_hgt2"]]
    a = store.annotations
    out = {
        "contig_taxa": sorted(c[rank].dropna().unique()),
        "region_taxa": sorted(R.loc[R["hgt1_rank"].notna(), f"r_{rank}"].dropna().unique()),
        "partner_taxa": sorted(E[f"p_{rank}"].dropna().unique()),
        "groups": sorted(int(g) for g in R["group_id"].dropna().unique()),
        "card_genes": [], "drug_classes": [], "rip_families": [],
    }
    if a is not None:
        card = a[a["source"] == "card"]
        out["card_genes"] = sorted(card["label"].dropna().unique())
        out["drug_classes"] = sorted({x for s in card["category"].dropna() for x in s.split(";") if x})
        out["rip_families"] = sorted(a.loc[a["source"] == "ripdb", "label"].dropna().unique())
    return out
