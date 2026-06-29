import streamlit as st
from data import dataset_sidebar

st.set_page_config(
    page_title="HGT Visualizer",
    layout="wide",
)

data_dir, assembly = dataset_sidebar()

st.title("HGT Analysis Visualizer")
st.markdown(
    f"""
**Assembly:** `{assembly or '(not set)'}`
**Data directory:** `{data_dir}`

This dashboard visualizes horizontal gene transfer (HGT) events detected using two complementary methods:

| Method | Files | Description |
|--------|-------|-------------|
| **HGT1** | `[d/p/c/o/f/g/s]HGT1.tsv` | Internal regions of contigs with taxonomy inconsistent with the contig background, inferred from taxonomic profiling of ORFs. The prefix indicates the taxonomic level at which the mismatch was detected. |
| **HGT1 DR** | `[...]HGT1_DR.tsv` | Donor–recipient contig pairs for Method 1 events where another contig exists with a taxonomy matching the transferred region. |
| **HGT2** | `[d/p/c/o/f/g/s]HGT2.tsv` | Identical open reading frames shared between contigs that have different taxonomic assignments. |

Navigate using the sidebar. Set `HGT_DATA_DIR` and `HGT_ASSEMBLY` environment variables to pre-configure the dataset at launch.
"""
)
