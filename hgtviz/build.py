#!/usr/bin/env python3
"""Build the hgt_viz store from a kairos-lr (e2e) results directory.

    python -m hgtviz.build --results path/to/results --out path/to/store \
        [--catalogue NAME] [--card-index aro_index.tsv]

Reads (under results/detect/, prefix = catalogue name):
  required  NAME.kairos_region_metadata.tsv
  optional  NAME.kairos_region_overlaps.tsv          shared ORFs per region pair
            NAME.assignRes2.tsv                      contig taxonomy
            NAME.potential_hgt_regions.taxonomy.tsv  region taxonomy
            NAME.hgt/NAME.{d..s}HGT1.tsv, allHGT2.tsv, allDR.tsv,
                     hgt_groups.tsv, hgt_report.tsv   get_hgt-v2.py calls
            NAME.features/events.tsv                 refined shared segments
            NAME.ffn_coords.tsv                      prodigal ORF coordinates
            NAME-card.tsv, NAME-ripdb.tsv            diamond best hits
            NAME.faa.eggs.emapper.annotations        eggNOG-mapper (2.0 or 2.1) output
  and results/call/insitu_events.tsv (optional).

Writes Parquet tables (contigs, regions, links, groups, dr, genes, annotations,
insitu, eggnog, cog_counts) plus meta.json into --out. eggnog.parquet is sorted by
contig and read per contig by the app, never loaded whole. The base table is `contigs`: every contig with
an HGT1 region or an HGT2 link.

Taxonomy is parsed here from the mmseqs tables rather than taken from the
get_hgt-v2.py outputs, whose lineage parsing mangles names such as 'Pseudomonas_A'.
The .hgt/ tables are used for membership only (which regions/pairs were called).
"""
import argparse
import datetime as dt
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

from . import RANKS
from . import cog as COG
from .taxonomy import coarsest, discordant_rank, read_mmseqs_taxonomy

STORE_VERSION = 1
DIAMOND_COLS = ["qtitle", "stitle", "pident", "bitscore", "evalue"]


# ── helpers ──────────────────────────────────────────────────────────────────

def _read_tsv(path, **kw):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except pd.errors.EmptyDataError:
        return None


def _join(values):
    vals = sorted({str(v) for v in values if pd.notna(v) and str(v) != ""})
    return ",".join(vals)


def _pair_key(a, b):
    a, b = np.asarray(a, dtype=object), np.asarray(b, dtype=object)
    lo = np.where(a <= b, a, b)
    hi = np.where(a <= b, b, a)
    return pd.Series(lo) + "\t" + pd.Series(hi)


def detect_catalogue(detect_dir):
    hits = sorted(glob.glob(os.path.join(detect_dir, "*.kairos_region_metadata.tsv")))
    if not hits:
        return None
    if len(hits) > 1:
        names = ", ".join(os.path.basename(h) for h in hits)
        raise SystemExit(f"Several catalogues in {detect_dir} ({names}); pass --catalogue")
    return os.path.basename(hits[0]).replace(".kairos_region_metadata.tsv", "")


# ── layers ───────────────────────────────────────────────────────────────────

def load_regions(det, log):
    regions = _read_tsv(f"{det}.kairos_region_metadata.tsv")
    if regions is None:
        raise SystemExit(f"Missing or empty {det}.kairos_region_metadata.tsv")
    regions = regions.rename(columns={"stop": "end"})
    keep = ["region_idx", "contig", "start", "end", "region_len", "contig_length",
            "region_type", "region_title", "orf_cts", "orfs_encompassed"]
    regions = regions[[c for c in keep if c in regions.columns]].copy()
    regions["region_type"] = regions["region_type"].str.strip("_")
    log(f"regions: {len(regions):,}")
    return regions


def load_region_taxonomy(det, regions, log):
    path = f"{det}.potential_hgt_regions.taxonomy.tsv"
    if not os.path.exists(path):
        log("region taxonomy: not found")
        return None
    tax = read_mmseqs_taxonomy(
        path, ids=regions["region_title"],
        id_transform=lambda s: s.str.replace(r"::.*$", "", regex=True))
    log(f"region taxonomy: {len(tax):,} regions")
    return tax


def load_contig_taxonomy(det, contigs, log):
    path = f"{det}.assignRes2.tsv"
    if not os.path.exists(path):
        log("contig taxonomy: not found")
        return None
    tax = read_mmseqs_taxonomy(path, ids=contigs)
    log(f"contig taxonomy: {len(tax):,} contigs")
    return tax


def hgt1_calls(hgt_dir, cat, regions, rtax, ctax, log):
    """region_idx -> rank of the HGT1 call. Taken from get_hgt-v2.py's per-rank
    tables when present, otherwise recomputed from the two taxonomies."""
    calls = {}
    if hgt_dir and os.path.isdir(hgt_dir):
        for r in RANKS:
            df = _read_tsv(os.path.join(hgt_dir, f"{cat}.{r}HGT1.tsv"), usecols=["region_idx"])
            if df is not None:
                for idx in df["region_idx"]:
                    calls.setdefault(idx, r)
        if calls:
            log(f"HGT1: {len(calls):,} regions (from {cat}.hgt/)")
            return pd.Series(calls, name="hgt1_rank"), "pipeline"
    if rtax is None or ctax is None:
        log("HGT1: no calls and no taxonomy to compute them")
        return pd.Series(dtype=object, name="hgt1_rank"), "none"
    m = (regions[["region_idx", "region_title", "contig"]]
         .merge(rtax[["id", "taxid"] + RANKS], left_on="region_title", right_on="id")
         .merge(ctax[["id", "taxid"] + RANKS], left_on="contig", right_on="id",
                suffixes=("_r", "_c")))
    m = m[(m["taxid_r"] > 1) & (m["taxid_c"] > 1)]
    rank = discordant_rank(m, m, "{}_r", "{}_c")
    s = pd.Series(rank.values, index=m["region_idx"].values, name="hgt1_rank").dropna()
    log(f"HGT1: {len(s):,} regions (recomputed from taxonomy)")
    return s, "computed"


def load_links(det, hgt_dir, cat, regions, log):
    """One row per region pair that shares sequence (undirected)."""
    ov = _read_tsv(f"{det}.kairos_region_overlaps.tsv")
    ev = _read_tsv(os.path.join(f"{det}.features", "events.tsv"))
    if ov is None and ev is None:
        log("links: no region overlaps or events")
        return pd.DataFrame(columns=["region_a", "region_b", "contig_a", "contig_b"])

    if ev is not None and len(ev):
        links = pd.DataFrame({
            "event_id": ev["event_id"],
            "region_a": ev["region_idx_1"], "contig_a": ev["contig_1"],
            "start_a": ev["shared_start_1"], "end_a": ev["shared_end_1"],
            "region_b": ev["region_idx_2"], "contig_b": ev["contig_2"],
            "start_b": ev["shared_start_2"], "end_b": ev["shared_end_2"],
            "strand": ev["strand"], "identity": ev["identity"],
            "shared_len": ev["shared_len"], "refined": ev["refined"],
        })
        src = "events"
    else:
        pairs = ov[["region_idx_contig1", "contig1", "region_idx_contig2", "contig2"]].drop_duplicates()
        pairs["key"] = _pair_key(pairs["region_idx_contig1"], pairs["region_idx_contig2"]).values
        pairs = pairs.drop_duplicates("key")
        links = pd.DataFrame({
            "region_a": pairs["region_idx_contig1"].values, "contig_a": pairs["contig1"].values,
            "region_b": pairs["region_idx_contig2"].values, "contig_b": pairs["contig2"].values,
        })
        reg = regions.set_index("region_idx")
        for s in "ab":
            links[f"start_{s}"] = links[f"region_{s}"].map(reg["start"])
            links[f"end_{s}"] = links[f"region_{s}"].map(reg["end"])
        src = "region overlaps"
    links = links[links["region_a"] != links["region_b"]].reset_index(drop=True)
    links["key"] = _pair_key(links["region_a"], links["region_b"]).values

    if ov is not None:
        ov["key"] = _pair_key(ov["region_idx_contig1"], ov["region_idx_contig2"]).values
        n_orfs = ov.drop_duplicates(["key", "orf_region1", "orf_region2"]).groupby("key").size()
        links["n_shared_orfs"] = links["key"].map(n_orfs).fillna(0).astype(int)

    # get_hgt-v2.py criteria-2 membership
    hgt2 = _read_tsv(os.path.join(hgt_dir, f"{cat}.allHGT2.tsv"),
                     usecols=["region_idx_contig1", "region_idx_contig2"]) if hgt_dir else None
    if hgt2 is not None:
        keys = set(_pair_key(hgt2["region_idx_contig1"], hgt2["region_idx_contig2"]))
        links["hgt2_called"] = links["key"].isin(keys)
    log(f"links: {len(links):,} region pairs (from {src})")
    return links.drop(columns="key")


def load_groups(hgt_dir, cat, log):
    if not hgt_dir:
        return None, None
    members = _read_tsv(os.path.join(hgt_dir, f"{cat}.hgt_groups.tsv"),
                        usecols=["region_idx", "hgt_group_id", "hgt_group_size"])
    report = _read_tsv(os.path.join(hgt_dir, f"{cat}.hgt_report.tsv"))
    if members is not None:
        log(f"groups: {members['hgt_group_id'].nunique():,} groups, {len(members):,} regions")
    return members, report


def load_dr(hgt_dir, cat, log):
    dr = _read_tsv(os.path.join(hgt_dir, f"{cat}.allDR.tsv")) if hgt_dir else None
    if dr is None or dr.empty:
        return pd.DataFrame(columns=["recipient", "recipient_region", "donor", "donor_region", "rank"])
    rank = pd.Series(np.nan, index=dr.index, dtype=object)
    for r in RANKS:
        col = f"{r}_shared_taxon"
        if col in dr.columns:
            rank = rank.where(rank.notna() | dr[col].isna(), r)
    out = pd.DataFrame({
        "recipient": dr["recipient"], "recipient_region": dr["recipient_region_idx"],
        "donor": dr["donor"], "donor_region": dr["donor_region_idx"], "rank": rank,
    }).drop_duplicates().reset_index(drop=True)
    log(f"donor-recipient pairs: {len(out):,}")
    return out


def load_genes(det, keep_contigs, regions, log):
    path = f"{det}.ffn_coords.tsv"
    df = _read_tsv(path, header=None, usecols=[0, 1, 2, 3],
                   names=["orf_id", "start", "end", "strand"], dtype={"orf_id": str})
    if df is None:
        log("genes: no ffn_coords.tsv")
        return None
    df["contig"] = df["orf_id"].str.rsplit("_", n=1).str[0]
    df = df[df["contig"].isin(keep_contigs)].copy()
    df["start"] = df["start"].astype(int) - 1   # prodigal is 1-based inclusive
    df["end"] = df["end"].astype(int)
    df["strand"] = pd.to_numeric(df["strand"], errors="coerce").fillna(0).astype(int)
    orf2region = orf_region_map(regions)
    df["region_idx"] = df["orf_id"].map(orf2region)
    log(f"genes: {len(df):,} on {df['contig'].nunique():,} contigs")
    return df.reset_index(drop=True)


def orf_region_map(regions):
    if "orfs_encompassed" not in regions.columns:
        return pd.Series(dtype=object)
    x = regions[["region_idx", "orfs_encompassed"]].dropna()
    x = x.assign(orf_id=x["orfs_encompassed"].str.split(",")).explode("orf_id")
    return x.drop_duplicates("orf_id").set_index("orf_id")["region_idx"]


def _read_diamond(path):
    df = _read_tsv(path, header=None, names=DIAMOND_COLS, dtype={"qtitle": str, "stitle": str},
                   quoting=3)
    if df is None:
        return None
    df["orf_id"] = df["qtitle"].str.split(" ", n=1).str[0]
    df["contig"] = df["orf_id"].str.rsplit("_", n=1).str[0]
    return df


def read_card_index(path):
    idx = _read_tsv(path, dtype=str)
    if idx is None:
        return None
    idx = idx.rename(columns=lambda c: c.strip())
    cols = {"ARO Accession": "accession", "AMR Gene Family": "family",
            "Drug Class": "category", "Resistance Mechanism": "mechanism"}
    missing = [c for c in cols if c not in idx.columns]
    if missing:
        raise SystemExit(f"{path}: missing columns {missing}")
    return idx[list(cols)].rename(columns=cols).drop_duplicates("accession").set_index("accession")


def load_annotations(det, keep_contigs, regions, card_index, log):
    frames = []
    card = _read_diamond(f"{det}-card.tsv")
    if card is not None:
        card = card[card["contig"].isin(keep_contigs)].copy()
        parts = card["stitle"].str.split("|")
        card["accession"] = parts.apply(lambda p: next((x for x in p if x.startswith("ARO:")), np.nan))
        card["label"] = parts.str[-1].str.replace(r"\s*\[.*\]\s*$", "", regex=True).str.strip()
        card["source"] = "card"
        for c in ["family", "category", "mechanism"]:
            card[c] = np.nan
        if card_index is not None:
            for c in ["family", "category", "mechanism"]:
                card[c] = card["accession"].map(card_index[c])
        frames.append(card)
        log(f"CARD hits on kept contigs: {len(card):,}"
            + ("" if card_index is not None else " (no aro_index: no drug classes)"))
    rip = _read_diamond(f"{det}-ripdb.tsv")
    if rip is not None:
        rip = rip[rip["contig"].isin(keep_contigs)].copy()
        fam = rip["stitle"].str.split(" ", n=1).str[0].str.rsplit("|", n=1)
        rip["label"] = fam.apply(lambda p: p[-1] if len(p) > 1 else "unknown")
        rip["accession"] = fam.str[0]
        rip["family"] = rip["label"]
        rip["category"] = np.where(rip["label"].isin(["parA", "parB"]), "partition", "replication")
        rip["mechanism"] = np.nan
        rip["source"] = "ripdb"
        frames.append(rip)
        log(f"RIP-DB hits on kept contigs: {len(rip):,}")
    if not frames:
        return None
    ann = pd.concat(frames, ignore_index=True)
    ann["region_idx"] = ann["orf_id"].map(orf_region_map(regions))
    ann["in_region"] = ann["region_idx"].notna()
    return ann[["orf_id", "contig", "source", "label", "accession", "family", "category",
                "mechanism", "pident", "bitscore", "evalue", "stitle", "region_idx", "in_region"]]


# store column: (emapper 2.1 name, emapper 2.0 name)
EGGNOG_COLS = {
    "evalue": ("evalue", "seed_ortholog_evalue"),
    "preferred_name": ("Preferred_name", "Preferred_name"),
    "description": ("Description", "eggNOG free text desc."),
    "cog_category": ("COG_category", "COG Functional cat."),
    "kegg_ko": ("KEGG_ko", "KEGG_ko"),
    "pfams": ("PFAMs", None),
    "cazy": ("CAZy", "CAZy"),
    "ec": ("EC", "EC"),
    "max_annot_lvl": ("max_annot_lvl", "best_tax_level"),
}


def find_eggnog(det):
    for path in [f"{det}.faa.eggs.emapper.annotations", f"{det}.emapper.annotations"]:
        if os.path.exists(path):
            return path
    hits = sorted(glob.glob(f"{det}*.emapper.annotations"))
    return hits[0] if hits else None


def _emapper_header(path):
    """Number of comment lines before the header, and the header's columns. emapper
    2.1 writes '##' comments and a '#query' header; 2.0 writes '#' comments and a
    '#query_name' header."""
    with open(path) as fh:
        for i, line in enumerate(fh):
            if line.startswith("#query"):
                return i, line[1:].rstrip("\n").split("\t")
            if not line.startswith("#"):
                break
    raise SystemExit(f"{path}: no '#query' header line; is this an emapper .annotations file?")


def load_eggnog(path, genes, log, chunksize=2_000_000):
    """Stream an emapper .annotations file, keeping only ORFs in the gene table and
    the columns in EGGNOG_COLS ('-' becomes missing)."""
    skip, header = _emapper_header(path)
    query = header[0]
    pick = {}
    for col, names in EGGNOG_COLS.items():
        name = next((n for n in names if n and n in header), None)
        if name:
            pick[name] = col
    keep = pd.Index(genes["orf_id"])
    parts, n_read = [], 0
    reader = pd.read_csv(path, sep="\t", skiprows=skip + 1, header=None, names=header,
                         usecols=[query] + list(pick), dtype=str, quoting=3,
                         na_values=["-"], keep_default_na=False, chunksize=chunksize)
    for chunk in reader:
        n_read += len(chunk)
        chunk = chunk[chunk[query].isin(keep)]   # also drops trailing '#' comment lines
        parts.append(chunk)
        log(f"  eggNOG: read {n_read:,} rows, kept {sum(len(x) for x in parts):,}")
    if not sum(len(x) for x in parts):
        log(f"WARNING: none of the {n_read:,} eggNOG queries match an ORF ID in ffn_coords.tsv "
            f"(e.g. {genes['orf_id'].iloc[0]!r}); was emapper run on a different .faa?")
    eg = pd.concat(parts, ignore_index=True).rename(columns={query: "orf_id", **pick})
    eg = eg.drop_duplicates("orf_id")
    for col in EGGNOG_COLS:
        if col not in eg.columns:
            eg[col] = np.nan
    eg["evalue"] = pd.to_numeric(eg["evalue"], errors="coerce")
    g = genes.set_index("orf_id")
    eg["contig"] = eg["orf_id"].map(g["contig"])
    eg["region_idx"] = eg["orf_id"].map(g["region_idx"])
    eg["in_region"] = eg["region_idx"].notna()
    eg["cog_group"] = eg["cog_category"].map(COG.group)
    eg = eg.sort_values(["contig", "orf_id"], kind="stable").reset_index(drop=True)
    log(f"eggNOG: {len(eg):,} of {len(genes):,} genes annotated "
        f"({eg['cog_category'].notna().sum():,} with a COG category)")
    return eg[["orf_id", "contig", "region_idx", "in_region", "preferred_name", "description",
               "cog_category", "cog_group", "kegg_ko", "pfams", "cazy", "ec", "evalue",
               "max_annot_lvl"]]


def cog_counts(genes, eg):
    """Genes per (contig, COG letter, inside a kairos region). Genes with several
    letters count once per letter; genes without a COG count under '-'."""
    g = genes[["orf_id", "contig", "region_idx"]].copy()
    g["cog"] = g["orf_id"].map(eg.set_index("orf_id")["cog_category"])
    g["cog"] = g["cog"].map(lambda c: COG.letters(c) or ["-"])
    g = g.explode("cog")
    g["in_region"] = g["region_idx"].notna()
    out = g.groupby(["contig", "cog", "in_region"], observed=True).size().rename("n").reset_index()
    return out


def load_insitu(results, log):
    df = _read_tsv(os.path.join(results, "call", "insitu_events.tsv"))
    if df is None:
        return None
    log(f"in situ events: {len(df):,}")
    return df


# ── main build ───────────────────────────────────────────────────────────────

def build(results, out, catalogue=None, card_index=None, eggnog=None, log=print):
    detect_dir = os.path.join(results, "detect")
    if not os.path.isdir(detect_dir):
        detect_dir = results  # allow pointing straight at results/detect
    cat = catalogue or detect_catalogue(detect_dir)
    if not cat:
        raise SystemExit(f"No *.kairos_region_metadata.tsv in {detect_dir}")
    det = os.path.join(detect_dir, cat)
    hgt_dir = f"{det}.hgt"
    if not os.path.isdir(hgt_dir):
        hgt_dir = None
        log("note: no .hgt/ directory (get_hgt-v2.py outputs); HGT calls recomputed from taxonomy")
    log(f"catalogue: {cat}")

    regions = load_regions(det, log)
    rtax = load_region_taxonomy(det, regions, log)
    links = load_links(det, hgt_dir, cat, regions, log)

    linked = set(links["contig_a"]) | set(links["contig_b"])
    ctax = load_contig_taxonomy(det, set(regions["contig"]) | linked, log)

    # region taxonomy columns r_<rank>
    if rtax is not None:
        rt = rtax.set_index("id")
        regions["region_taxon"] = regions["region_title"].map(rt["tax_name"])
        for r in RANKS:
            regions[f"r_{r}"] = regions["region_title"].map(rt[r])
    else:
        regions["region_taxon"] = np.nan
        for r in RANKS:
            regions[f"r_{r}"] = np.nan

    hgt1, hgt1_source = hgt1_calls(hgt_dir, cat, regions, rtax, ctax, log)
    regions["hgt1_rank"] = regions["region_idx"].map(hgt1)

    # contig taxonomy on both link ends; HGT2 = discordant contig lineages
    if ctax is not None:
        ct = ctax.set_index("id")
        for s in "ab":
            for r in RANKS:
                links[f"{s}_{r}"] = links[f"contig_{s}"].map(ct[r])
        links["discordant_rank"] = discordant_rank(links, links, "a_{}", "b_{}")
    else:
        for s in "ab":
            for r in RANKS:
                links[f"{s}_{r}"] = np.nan
        links["discordant_rank"] = np.nan
    if "hgt2_called" in links.columns:
        links["is_hgt2"] = links["hgt2_called"]
        hgt2_source = "pipeline"
    else:
        links["is_hgt2"] = links["discordant_rank"].notna()
        hgt2_source = "computed" if ctax is not None else "none"

    members, report = load_groups(hgt_dir, cat, log)
    if members is not None:
        regions = regions.merge(members.rename(columns={"hgt_group_id": "group_id",
                                                        "hgt_group_size": "group_size"}),
                                on="region_idx", how="left")
    else:
        regions["group_id"] = np.nan
        regions["group_size"] = np.nan

    dr = load_dr(hgt_dir, cat, log)

    # ── base: contigs with HGT1 or HGT2 ──────────────────────────────────────
    hgt1_contigs = set(regions.loc[regions["hgt1_rank"].notna(), "contig"])
    h2 = links[links["is_hgt2"]]
    hgt2_contigs = set(h2["contig_a"]) | set(h2["contig_b"])
    base = sorted(hgt1_contigs | hgt2_contigs)
    log(f"base contigs: {len(base):,} (HGT1 {len(hgt1_contigs):,}, HGT2 {len(hgt2_contigs):,})")

    # keep every link touching a base contig, and regions/genes on both ends
    links = links[links["contig_a"].isin(base) | links["contig_b"].isin(base)].reset_index(drop=True)
    gene_contigs = set(base) | set(links["contig_a"]) | set(links["contig_b"])
    regions = regions[regions["contig"].isin(gene_contigs)].reset_index(drop=True)

    contigs = pd.DataFrame({"contig": base})
    length = regions.groupby("contig")["contig_length"].first()
    contigs["length"] = contigs["contig"].map(length)
    if ctax is not None:
        ct = ctax.set_index("id")
        contigs["taxon"] = contigs["contig"].map(ct["tax_name"])
        contigs["tax_level"] = contigs["contig"].map(ct["tax_level"])
        contigs["tax_agreement"] = contigs["contig"].map(ct["percent_agreed"])
        for r in RANKS:
            contigs[r] = contigs["contig"].map(ct[r])
    else:
        for c in ["taxon", "tax_level", "tax_agreement"] + RANKS:
            contigs[c] = np.nan

    rg = regions.groupby("contig")
    contigs["n_regions"] = contigs["contig"].map(rg.size()).fillna(0).astype(int)
    h1 = regions[regions["hgt1_rank"].notna()].groupby("contig")
    contigs["hgt1"] = contigs["contig"].isin(hgt1_contigs)
    contigs["n_hgt1_regions"] = contigs["contig"].map(h1.size()).fillna(0).astype(int)
    contigs["hgt1_rank"] = contigs["contig"].map(h1["hgt1_rank"].agg(coarsest))

    ends = pd.concat([
        h2[["contig_a", "contig_b", "discordant_rank"]].set_axis(["contig", "partner", "rank"], axis=1),
        h2[["contig_b", "contig_a", "discordant_rank"]].set_axis(["contig", "partner", "rank"], axis=1),
    ])
    eg = ends.groupby("contig")
    contigs["hgt2"] = contigs["contig"].isin(hgt2_contigs)
    contigs["n_hgt2_partners"] = contigs["contig"].map(eg["partner"].nunique()).fillna(0).astype(int)
    contigs["hgt2_rank"] = contigs["contig"].map(eg["rank"].agg(coarsest))

    contigs["groups"] = contigs["contig"].map(rg["group_id"].agg(
        lambda s: _join(s.dropna().astype(int)))).fillna("")
    contigs["n_groups"] = contigs["groups"].str.split(",").apply(lambda x: len([y for y in x if y]))

    rec, don = set(dr["recipient"]), set(dr["donor"])
    contigs["dr_role"] = [
        "both" if c in rec and c in don else "recipient" if c in rec else "donor" if c in don else ""
        for c in contigs["contig"]]

    genes = load_genes(det, gene_contigs, regions, log)
    ann = load_annotations(det, gene_contigs,
                           pd.concat([regions], ignore_index=True),
                           read_card_index(card_index) if card_index else None, log)
    for src, col in [("card", "card_genes"), ("ripdb", "rip_families")]:
        if ann is not None and (ann["source"] == src).any():
            a = ann[ann["source"] == src]
            contigs[col] = contigs["contig"].map(a.groupby("contig")["label"].agg(_join)).fillna("")
            contigs[f"{col}_in_region"] = contigs["contig"].map(
                a[a["in_region"]].groupby("contig")["label"].agg(_join)).fillna("")
        else:
            contigs[col] = ""
            contigs[f"{col}_in_region"] = ""
    if ann is not None and "category" in ann:
        a = ann[(ann["source"] == "card") & ann["category"].notna()]
        contigs["drug_classes"] = contigs["contig"].map(
            a.groupby("contig")["category"].agg(lambda s: _join(";".join(s).split(";")))).fillna("")
    else:
        contigs["drug_classes"] = ""

    eggnog_path = eggnog or find_eggnog(det)
    eg = cc = None
    if eggnog_path and genes is not None:
        log(f"eggNOG: {eggnog_path}")
        eg = load_eggnog(eggnog_path, genes, log)
        cc = cog_counts(genes[genes["contig"].isin(set(contigs["contig"]))], eg)  # base contigs only
        e = eg[eg["cog_category"].notna()]
        contigs["cogs"] = contigs["contig"].map(
            e.groupby("contig")["cog_category"].agg(lambda s: "".join(sorted(set("".join(s)))))).fillna("")
    elif eggnog_path:
        log("eggNOG: skipped, needs the gene table (ffn_coords.tsv)")

    insitu = load_insitu(results, log)
    if insitu is not None and {"recipient_contig", "donor_contig"} <= set(insitu.columns):
        contigs["insitu_as_recipient"] = contigs["contig"].map(
            insitu["recipient_contig"].value_counts()).fillna(0).astype(int)
        contigs["insitu_as_donor"] = contigs["contig"].map(
            insitu["donor_contig"].value_counts()).fillna(0).astype(int)

    # ── write ────────────────────────────────────────────────────────────────
    os.makedirs(out, exist_ok=True)
    tables = {"contigs": contigs, "regions": regions.drop(columns=["orfs_encompassed"], errors="ignore"),
              "links": links, "dr": dr, "groups": report, "genes": genes,
              "annotations": ann, "insitu": insitu, "eggnog": eg, "cog_counts": cc}
    written = {}
    for name, df in tables.items():
        path = os.path.join(out, f"{name}.parquet")
        if df is None:
            if os.path.exists(path):
                os.remove(path)
            continue
        if name == "eggnog":   # small row groups: the app reads one contig at a time
            df.to_parquet(path, index=False, row_group_size=20_000)
        else:
            df.to_parquet(path, index=False)
        written[name] = len(df)
    meta = {
        "store_version": STORE_VERSION,
        "catalogue": cat,
        "results": os.path.abspath(results),
        "built": dt.datetime.now().isoformat(timespec="seconds"),
        "tables": written,
        "hgt1_source": hgt1_source,
        "hgt2_source": hgt2_source,
        "card_index": os.path.abspath(card_index) if card_index else None,
        "eggnog": os.path.abspath(eggnog_path) if eg is not None else None,
    }
    with open(os.path.join(out, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    log(f"wrote {out}: " + ", ".join(f"{k} {v:,}" for k, v in written.items()))
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="kairos-lr results/ directory (or results/detect)")
    ap.add_argument("--out", required=True, help="output store directory")
    ap.add_argument("--catalogue", help="catalogue name (auto-detected if only one)")
    ap.add_argument("--card-index", help="CARD aro_index.tsv, for drug class / gene family / mechanism")
    ap.add_argument("--eggnog", help="emapper .annotations file (default: NAME.faa.eggs.emapper.annotations "
                                     "next to the other detect outputs)")
    a = ap.parse_args(argv)
    build(a.results, a.out, a.catalogue, a.card_index, a.eggnog)


if __name__ == "__main__":
    sys.exit(main())
