# HGT viz

Interactive Streamlit browser for horizontal gene transfer (HGT) calls from the
[kairos-lr](https://github.com/clb21565/kairos-lr) `e2e` pipeline.

The unit is the **contig**: every contig with an HGT1 region (region taxonomy differs
from the contig) or an HGT2 link (a shared region with a contig of a different
lineage). You filter that list, then select a contig to see the layers on top of it:
its regions, HGT2 partners, HGT groups, donor–recipient pairs, CARD / RIP-DB hits and
in situ calls.

## For colleagues: try the app

The app runs on a compute node of the HPC and you view it in your own browser through
an SSH tunnel. Nothing needs installing on your computer.

1. Log in to the HPC and run the launcher (the full path is in the message you got
   with this link):

   ```bash
   bash /abs/path/to/hgt_viz/slurm/launch.sh
   ```

   It submits a small job (2 CPUs, 32 GB, 8 h), waits for it to start and prints an
   `ssh -N -L …` command. If it says port 8501 is busy, rerun with `PORT=8642 bash …`.
2. On **your own computer**, paste that `ssh` command into a terminal and leave it
   open (Windows: PowerShell works).
3. Open <http://localhost:8501> (or the port you chose). The first page takes about a
   minute while the data loads.
4. When you are done: `scancel <jobid>` (also printed by the launcher).

### A five-minute tour

- **HGT calls explained** (left menu): what the two call types mean. Start here.
- **Contig browser**: pick a contig in the table; the gene map, the lineage table and
  the tabs below describe it. Try the *Show contigs with* switch above the table and
  the *Evidence shown* switch above the map.
- **Function search**: e.g. eggNOG text `transposase` or a CARD gene, *Where* =
  *inside a region with an HGT call*, **Apply**; select a gene and open its contig
  (the gene is marked in red).
- **Sidebar → Taxonomy → Search full lineages**: e.g. host `Burkholderiaceae`,
  foreign/partner `Comamonas`.
- **Summary** and **HGT groups**: overview of whatever the sidebar filters select.

### Feedback we would like

- Which calls look implausible, and what in the app made you think so (contig + region
  IDs help)?
- What did you look for and not find? Which labels or colours were unclear?
- Anything slow (which page, which action)?

### Setting it up for others (maintainer)

Copy `slurm/site.env.example` to `slurm/site.env` (ignored by git) and fill in absolute
paths: the store, the results it came from, the conda env (a path colleagues can read,
e.g. `/home/<you>/.conda/envs/hgt_viz`) and the login host. Make the repo, store and env
readable for your colleagues (e.g. `chmod -R g+rX` or `o+rX` on those directories and
`o+x` on each parent directory). Then send them the absolute path of `slurm/launch.sh`.

## Pages

| Page | What it shows |
|------|---------------|
| **Contig browser** | Filterable contig table (quick switch: any / HGT1 / HGT2 / both / one type only), TSV exports of the filtered set, and for the selected contig: a gene map (regions, shared segments; genes coloured by HGT evidence, CARD / RIP-DB hit, eggNOG COG group or identical ORFs; CARD or eggNOG gene names under the genes where they fit; wrapped eggNOG text on hover; an **Evidence** switch to show HGT1 or HGT2 only; the gene opened from Function search is marked in red and the map zooms to its region), a **lineage table** aligning the host's lineage with each foreign region's and each HGT2 partner lineage (grouped by the diverging taxon, or full lineages), shaded where they first differ, tabs for regions, partners (select a row for an aligned two-contig view), donor–recipient pairs, groups, annotations (CARD / RIP-DB and eggNOG) and in situ calls |
| **Summary** | For the filtered contigs: rank of the HGT calls, HGT1 region → host taxon matrix, HGT2 contig ↔ partner matrix, COG categories inside vs. outside kairos regions, CARD genes / drug classes / RIP-DB families inside vs. outside kairos regions |
| **Function search** | Search genes by CARD gene, drug class, RIP-DB family, COG category or eggNOG text: anywhere on the contig, inside a kairos region, or inside a region with an HGT call (optionally at chosen ranks). Each matching gene shows its region's HGT level (coarsest rank), HGT1 rank, HGT2 rank and number of HGT2 partners; a count splits hits into HGT regions / other kairos regions / outside regions. Shows matching genes in the whole store and on the filtered contigs; **Apply as filter** keeps only contigs with a match on every page; select a gene to open its contig |
| **HGT groups** | Groups that touch the filtered contigs, their CARD / RIP-DB content, member regions, and a button to filter the browser to one group |
| **HGT calls explained** | What HGT1 (*foreign region*) and HGT2 (*cross-lineage share*) mean, how they combine with donor–recipient, groups and in situ calls, and how to read the gene-map colours |

Display names: HGT1 is shown as **Foreign region** (the region's taxonomy differs from its
host contig) and HGT2 as **Cross-lineage share** (the same region on contigs of different
lineages). Change them in `hgtviz/__init__.py`; pipeline files keep HGT1 / HGT2.

Sidebar filters apply to every page: contig ID, required evidence (HGT1, HGT2,
donor–recipient, in situ), host / HGT1-region / HGT2-partner taxon at any rank, free-text
search in the full lineage of the host and of its foreign regions / partners (commas =
alternatives, a trailing `;` = exact name, e.g. host `Burkholderiaceae`, partner
`Comamonas`), rank of the HGT1 and HGT2 calls, HGT group and group size, contig length,
and the function filter set on the Function search page.

## 1. Build a store

The app reads a small directory of Parquet tables built once from a kairos-lr
`results/` directory:

```bash
python -m hgtviz.build \
    --results /path/to/kairos-lr/results \
    --out     /path/to/hgt_viz_store \
    --card-index /path/to/card/aro_index.tsv     # optional: drug class, gene family, mechanism
```

Inputs (under `results/detect/`, or directly in `--results` if it has no `detect/`;
prefix = catalogue name, auto-detected):

| File | Required | Used for |
|------|----------|----------|
| `NAME.kairos_region_metadata.tsv` | yes | regions |
| `NAME.kairos_region_overlaps.tsv` | | shared ORFs per region pair |
| `NAME.features/events.tsv` | | refined shared segments (coordinates, identity, strand) |
| `NAME.assignRes2.tsv`, `NAME.potential_hgt_regions.taxonomy.tsv` | | contig and region taxonomy |
| `NAME.hgt/` (get_hgt-v2.py outputs) | | HGT1 / HGT2 calls, groups, donor–recipient pairs |
| `NAME.ffn_coords.tsv` | | gene coordinates |
| `NAME-card.tsv`, `NAME-ripdb.tsv` | | diamond best hits (`qtitle stitle pident bitscore evalue`) |
| `NAME.faa.eggs.emapper.annotations` | | eggNOG-mapper 2.0 / 2.1 output (`--eggnog` to point elsewhere) |
| `NAME.ffn.clusts_cluster.tsv` | | mmseqs ORF clusters used by kairos-lr (`--orf-clusters` to point elsewhere) |
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
- Gene maps open in **HGT evidence** colouring: a gene is HGT if its region has an
  HGT1 call or it overlaps an HGT2 shared segment; its colour is the coarsest rank of
  disagreement over both (one blue ramp: domain darkest, then phylum, class, order,
  family-or-below). Regions carry labels like `HGT1 phylum · HGT2 phylum (13)`; hover
  gives the taxa (region vs host, partner lineages), partner count and best identity.
  Shared-segment ticks use the same ramp.
- ORF clusters: two streaming passes give each store gene its cluster and the
  cluster's catalogue-wide size (copies) and spread (contigs). The gene map's
  "Identical ORFs" mode colours at most 4 clusters per view, numbered on the map (the
  four hues that stay distinguishable pairwise): in a single contig, the clusters with
  copies on the most other contigs; in the two-contig view, clusters on both contigs,
  joined by ribbons. Other genes with an identical copy elsewhere are dark grey,
  unique genes light grey; hover shows every gene's cluster.
- eggNOG: the file is streamed in chunks and reduced to the genes in the store and the
  columns name, description, COG category, KEGG KO, Pfam, CAZy, EC, e-value (GO terms
  and KEGG pathways are dropped). `eggnog.parquet` is sorted by contig and the app reads
  it one contig at a time; it is never loaded whole. `cog_counts.parquet` (genes per
  base contig × COG letter × inside/outside a region) drives the COG filter and the
  Summary chart.

## 2. Run the app

```bash
conda create -n hgt_viz python=3.11 -y && conda activate hgt_viz
pip install -r requirements.txt
HGT_VIZ_STORE=/path/to/hgt_viz_store streamlit run app.py --server.headless true --server.port 8501
```

The store path can also be changed in the sidebar. Pointing it at a `results/`
directory offers to build the store there.

### On an HPC compute node

`slurm/hgt_viz.sbatch` builds the store (if `STORE/meta.json` is missing, or with
`REBUILD=1`) and then serves the app. Create the env once with
`conda env create -f environment.yml`, then:

```bash
sbatch --export=ALL,REPO=$PWD,RESULTS=/path/to/kairos-lr/results,STORE=/path/to/hgt_viz_store,\
CARD_INDEX=/path/to/aro_index.tsv slurm/hgt_viz.sbatch
```

Always pass the same `STORE` (it defaults to `hgt_viz_store/` next to `RESULTS`, which
is easy to lose track of). Add `REBUILD=1` after changes to the kairos-lr outputs or
after an hgt_viz update that adds store tables (the app shows which layers a store has:
e.g. no "Identical ORFs" colour option means the store predates ORF clusters). Build
memory and time go to `build_time_<jobid>.txt`. The log prints the `ssh -L` tunnel
command. `slurm/launch.sh` does the same for an existing store with settings from
`slurm/site.env` (see *For colleagues*).

After `git pull`, restart the job (`scancel`, then submit again); Streamlit does not
reload the `hgtviz` package while it runs.

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
