import os
import subprocess
import sys

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from hgtviz import filters as F  # noqa: E402
from hgtviz import store as S  # noqa: E402
from hgtviz.build import build  # noqa: E402
from hgtviz.taxonomy import discordant_rank, parse_lineage  # noqa: E402


def test_parse_lineage_keeps_gtdb_suffixes():
    lin = "d_Bacteria;p_Bacillota_A;c_Clostridia;g_Pseudomonas_A;s_Pseudomonas_A sp1"
    t = parse_lineage(lin)
    assert t["p"] == "Bacillota_A"
    assert t["g"] == "Pseudomonas_A"
    assert t["s"] == "Pseudomonas_A sp1"


def test_parse_lineage_skips_root_and_empty():
    assert parse_lineage("-_root;d_Archaea") == {"d": "Archaea"}
    assert parse_lineage("") == {}
    assert parse_lineage(None) == {}


def test_discordant_rank():
    a = pd.DataFrame({"d": ["B", "B", "B", "B"], "p": ["X", "X", "X", None],
                      "c": ["c1", "c1", "c2", None], "o": [None] * 4, "f": [None] * 4,
                      "g": [None] * 4, "s": [None] * 4})
    b = pd.DataFrame({"d": ["A", "B", "B", "B"], "p": ["X", "Y", "X", "Y"],
                      "c": ["c1", "c1", "c1", "c9"], "o": [None] * 4, "f": [None] * 4,
                      "g": [None] * 4, "s": [None] * 4})
    r = discordant_rank(a, b).tolist()
    assert r[:3] == ["d", "p", "c"]
    assert pd.isna(r[3])  # p unnamed on one side: never becomes discordant


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo")
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "make_demo.py"),
                    "--out", str(out), "--contigs", "200", "--elements", "25"],
                   check=True, stdout=subprocess.DEVNULL)
    store_dir = out / "store"
    build(str(out / "results"), str(store_dir), card_index=str(out / "aro_index.tsv"),
          log=lambda *_: None)
    return S.load(str(store_dir))


def test_build_tables(demo):
    c = demo.contigs
    assert len(c) > 0 and c["contig"].is_unique
    assert (c["hgt1"] | c["hgt2"]).all()          # base = HGT1 or HGT2
    assert demo.meta["hgt1_source"] == "computed"  # no .hgt/ dir in this demo
    # every HGT1 contig has an HGT1 region, every HGT2 contig an HGT2 link
    r = demo.regions
    assert set(c.loc[c["hgt1"], "contig"]) == set(r.loc[r["hgt1_rank"].notna(), "contig"])
    e = demo.link_ends
    assert set(c.loc[c["hgt2"], "contig"]) == set(e.loc[e["is_hgt2"], "contig"])
    # links are undirected: no region pair appears twice
    L = demo.links
    key = L.apply(lambda x: tuple(sorted((x["region_a"], x["region_b"]))), axis=1)
    assert key.is_unique


def test_annotations(demo):
    a = demo.annotations
    assert set(a["source"]) == {"card", "ripdb"}
    assert a.loc[a["source"] == "card", "category"].notna().all()   # aro_index joined
    assert a.loc[a["source"] == "card", "accession"].str.startswith("ARO:").all()
    assert not a.loc[a["source"] == "ripdb", "label"].eq("unknown").any()
    # in_region agrees with the gene table's region assignment
    g = demo.genes.set_index("orf_id")["region_idx"]
    assert (a["region_idx"].fillna("") == a["orf_id"].map(g).fillna("")).all()


def test_filters(demo):
    every = F.apply(demo, F.FilterSpec())
    assert len(every) == len(demo.contigs)
    only2 = F.apply(demo, F.FilterSpec(require=["HGT2"]))
    assert only2["hgt2"].all()
    card = demo.annotations[demo.annotations["source"] == "card"]
    gene = card["label"].value_counts().index[0]
    anywhere = F.apply(demo, F.FilterSpec(card_genes=[gene]))
    inside = F.apply(demo, F.FilterSpec(card_genes=[gene], function_scope="region"))
    assert set(inside["contig"]) <= set(anywhere["contig"])
    assert set(anywhere["contig"]) == set(card.loc[card["label"] == gene, "contig"]) & set(demo.contigs["contig"])
    opts = F.options(demo, "p")
    if opts["partner_taxa"]:
        t = opts["partner_taxa"][0]
        hit = F.apply(demo, F.FilterSpec(partner_taxa=[t], tax_rank="p"))
        e = demo.link_ends
        assert set(hit["contig"]) == set(e.loc[e["is_hgt2"] & (e["p_p"] == t), "contig"])


def test_subset_matches_isin(demo):
    """Coded subsets select the same rows as a string isin on contig names."""
    drug = sorted(F.options(demo, "p")["drug_classes"])[:1]
    for spec in [F.FilterSpec(), F.FilterSpec(require=["HGT2"]), F.FilterSpec(drug_classes=drug)]:
        filt = F.apply(demo, spec)
        names = set(filt["contig"])
        for name in ["regions", "link_ends", "annotations"]:
            df = demo.link_ends if name == "link_ends" else getattr(demo, name)
            sub = demo.subset(name, filt)
            assert sub.index.equals(df[df["contig"].isin(names)].index), (spec, name)
            assert (demo.contigs["contig"].to_numpy()[sub["_cid"]] == sub["contig"].to_numpy()).all()
        E = demo.subset("link_ends", filt, hgt2_only=True)
        assert (demo.ids()[E["_pid"]] == E["partner"].to_numpy()).all()
        assert E["is_hgt2"].all()
    # natural sort at load
    ids = demo.contigs["contig"].tolist()
    assert ids.index("contig_2") < ids.index("contig_10") if {"contig_2", "contig_10"} <= set(ids) else True


def test_eggnog(demo):
    from hgtviz import cog as COG
    assert demo.has("eggnog") and demo.cog_counts is not None
    c = demo.contigs["contig"].iloc[0]
    e = demo.eggnog_for(c)
    assert len(e) and (e["contig"] == c).all()
    # cog_counts: one row-sum per gene per letter, all genes accounted for
    cc = demo.cog_counts
    base = demo.genes[demo.genes["contig"].isin(set(demo.contigs["contig"]))]   # counts: base contigs
    per_gene = base.merge(pd.read_parquet(demo.eggnog_path, columns=["orf_id", "cog_category"]),
                                on="orf_id", how="left")
    n_letters = per_gene["cog_category"].map(lambda x: max(len(COG.letters(x)), 1)).sum()
    assert cc["n"].sum() == n_letters
    # COG filter and text search agree with the eggNOG table
    eg = pd.read_parquet(demo.eggnog_path)
    eg = eg[eg["contig"].isin(set(demo.contigs["contig"]))]   # eggnog also covers partners
    hit = F.apply(demo, F.FilterSpec(cog_letters=["X"]))
    assert set(hit["contig"]) == set(eg.loc[eg["cog_category"].str.contains("X", na=False), "contig"])
    t = F.apply(demo, F.FilterSpec(eggnog_text="TRANSPOSASE"))
    assert set(t["contig"]) == set(eg.loc[eg["description"].str.contains("ransposase", na=False), "contig"])
    t_in = F.apply(demo, F.FilterSpec(eggnog_text="transposase", function_scope="region"))
    assert set(t_in["contig"]) <= set(t["contig"])
    assert COG.group("KX") == COG.MOBILE and COG.group("SK") == COG.INFO
    assert COG.group(None) == COG.NONE and COG.group("RS") == COG.POOR


def test_emapper_v2_0_header(tmp_path):
    from hgtviz.build import load_eggnog
    path = tmp_path / "x.emapper.annotations"
    cols = ["query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue", "seed_ortholog_score",
            "best_tax_level", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway",
            "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
            "BiGG_Reaction", "taxonomic scope", "eggNOG OGs", "best eggNOG OG",
            "COG Functional cat.", "eggNOG free text desc."]
    row = ["c1_1", "x", "1e-30", "99", "Bacteria", "tnpA", "-", "-", "ko:K07483"] + ["-"] * 8 \
        + ["bact", "-", "-", "X", "Transposase"]
    path.write_text("# emapper 2.0\n# x\n#" + "\t".join(cols) + "\n" + "\t".join(row) + "\n"
                    + "# done\n")
    genes = pd.DataFrame({"orf_id": ["c1_1", "c1_2"], "contig": ["c1", "c1"],
                          "region_idx": ["r1", None]})
    eg = load_eggnog(str(path), genes, log=lambda *_: None)
    assert eg.loc[0, "preferred_name"] == "tnpA" and eg.loc[0, "cog_category"] == "X"
    assert eg.loc[0, "description"] == "Transposase" and bool(eg.loc[0, "in_region"])


def test_orf_clusters(demo):
    from hgtviz import genemap as G
    g = demo.genes
    assert {"cluster", "cl_members", "cl_contigs"} <= set(g.columns)
    assert (g["cluster"] >= 0).all() and (g["cl_contigs"] >= 1).all()
    assert (g["cl_members"] >= g["cl_contigs"]).all()
    # copies of an element gene on two linked contigs share a cluster
    e = demo.link_ends[demo.link_ends["is_hgt2"]].iloc[0]
    g1, g2 = demo.rows_for("genes", e["contig"]), demo.rows_for("genes", e["partner"])
    shared = set(g1["cluster"]) & set(g2["cluster"])
    assert shared
    a, b = G.gene_category(g1, None), G.gene_category(g2, None)
    pick = G.assign_clusters([a, b], paired=True)
    assert 0 < len(pick) <= len(G.CLUSTER_COLORS) and set(pick) <= shared
    assert set(a.loc[a["cluster"].isin(pick), "clcat"]) <= {f"cl{i}" for i in range(1, 5)}
    assert (a.loc[a["cl_contigs"] == 1, "clcat"] == "unique").all()


def test_missing_store(tmp_path):
    with pytest.raises(S.StoreError):
        S.load(str(tmp_path))
