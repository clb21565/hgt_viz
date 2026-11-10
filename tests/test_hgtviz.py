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


def test_missing_store(tmp_path):
    with pytest.raises(S.StoreError):
        S.load(str(tmp_path))
