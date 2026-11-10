# HGT viz

Interactive Streamlit browser for horizontal gene transfer (HGT) calls from the
[kairos-lr](https://github.com/clb21565/kairos-lr) `e2e` pipeline.

The unit is the **contig**: every contig with an HGT1 region (region taxonomy differs
from the contig) or an HGT2 link (a shared region with a contig of a different
lineage). You filter that list, then select a contig to see the layers on top of it:
its regions, HGT2 partners, HGT groups, donor–recipient pairs, CARD / RIP-DB hits and
in situ calls.

## Pages

| Page | What it shows |
|------|---------------|
| **Contig browser** | Filterable contig table, TSV exports of the filtered set, and for the selected contig: a gene map (regions, CARD and RIP-DB hits, shared segments), tabs for regions, partners (with an aligned two-contig view), donor–recipient pairs, groups, annotations and in situ calls |
| **Summary** | For the filtered contigs: rank of the HGT calls, HGT1 region → host taxon matrix, HGT2 contig ↔ partner matrix, CARD genes / drug classes / RIP-DB families inside vs. outside kairos regions |
| **HGT groups** | Groups that touch the filtered contigs, their CARD / RIP-DB content, member regions, and a button to filter the browser to one group |

Sidebar filters apply to every page: contig ID, required evidence (HGT1, HGT2,
donor–recipient, in situ), host / HGT1-region / HGT2-partner taxon at any rank, rank
of the HGT1 and HGT2 calls, HGT group and group size, contig length, CARD gene, CARD
drug class and RIP-DB family (anywhere on the contig, or only inside a kairos region).

## 1. Build a store

The app reads a small directory of Parquet tables built once from a kairos-lr
`results/` directory:

```bash
python -m hgtviz.build \
    --results /path/to/kairos-lr/results \
    --out     /path/to/hgt_viz_store \
    --card-index /path/to/card/aro_index.tsv     # optional: drug class, gene family, mechanism
```

Inputs (under `results/detect/`, prefix = catalogue name, auto-detected):

| File | Required | Used for |
|------|----------|----------|
| `NAME.kairos_region_metadata.tsv` | yes | regions |
| `NAME.kairos_region_overlaps.tsv` | | shared ORFs per region pair |
| `NAME.features/events.tsv` | | refined shared segments (coordinates, identity, strand) |
| `NAME.assignRes2.tsv`, `NAME.potential_hgt_regions.taxonomy.tsv` | | contig and region taxonomy |
| `NAME.hgt/` (get_hgt-v2.py outputs) | | HGT1 / HGT2 calls, groups, donor–recipient pairs |
| `NAME.ffn_coords.tsv` | | gene coordinates |
| `NAME-card.tsv`, `NAME-ripdb.tsv` | | diamond best hits (`qtitle stitle pident bitscore evalue`) |
| `results/call/insitu_events.tsv` | | in situ calls |

Notes:

- Without `NAME.hgt/`, HGT1 and HGT2 are recomputed from the two taxonomy tables with
  the same rule as `get_hgt-v2.py` (highest rank at which both lineages are named and
  differ). `meta.json` records which source was used.
- Taxonomy is parsed from the mmseqs tables directly: `get_hgt-v2.py` strips every
  `<letter>_` in a lineage, which turns GTDB names like `Pseudomonas_A` into
  `PseudomonaA`.
- RIP-DB family is the text after the last `|` in the subject header (RIP-DB v1.01).
  `parA` / `parB` hits are labelled `partition`, all others `replication`.
- Only genes and hits on base contigs and their link partners are kept, so the store
  stays small even for a large catalogue.

## 2. Run the app

```bash
conda create -n hgt_viz python=3.11 -y && conda activate hgt_viz
pip install -r requirements.txt
HGT_VIZ_STORE=/path/to/hgt_viz_store streamlit run app.py --server.headless true --server.port 8501
```

The store path can also be changed in the sidebar. Pointing it at a `results/`
directory offers to build the store there.

### On an HPC compute node

```bash
#!/bin/bash
#SBATCH --job-name=hgt_viz
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=hgt_viz_%j.log

conda activate hgt_viz
echo "Tunnel: ssh -L 8501:$(hostname):8501 ${USER}@<login-host>"
HGT_VIZ_STORE=/path/to/hgt_viz_store \
    streamlit run /path/to/hgt_viz/app.py --server.headless true --server.port 8501
```

Then run the printed `ssh -L` command locally and open http://localhost:8501.

## Demo data

`scripts/make_demo.py` simulates a small `results/` directory (synthetic taxa,
contigs, mobile elements, CARD and RIP-DB hits; the CARD accessions are made up):

```bash
python scripts/make_demo.py --out demo --get-hgt /path/to/kairos-lr/scripts/get_hgt-v2.py
python -m hgtviz.build --results demo/results --out demo/store --card-index demo/aro_index.tsv
HGT_VIZ_STORE=demo/store streamlit run app.py
```

Leave out `--get-hgt` to test the path where HGT calls are recomputed from taxonomy.

## Tests

```bash
python -m pytest tests
```
