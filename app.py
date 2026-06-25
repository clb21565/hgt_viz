import streamlit as st

st.set_page_config(
    page_title="HGT Visualizer",
    layout="wide",
)

st.title("HGT Analysis Visualizer")
st.markdown(
    """
**Assembly:** `damhusaen_as_rp3`

This dashboard visualizes horizontal gene transfer (HGT) events detected using two complementary methods:

| Method | Files | Description |
|--------|-------|-------------|
| **HGT1** | `[d/p/c/o/f/g/s]HGT1.tsv` | Internal regions of contigs with taxonomy inconsistent with the contig background, inferred from taxonomic profiling of ORFs. The prefix indicates the taxonomic level at which the mismatch was detected. |
| **HGT1 DR** | `[...]HGT1_DR.tsv` | Donor–recipient contig pairs for Method 1 events where another contig exists with a taxonomy matching the transferred region. |
| **HGT2** | `[d/p/c/o/f/g/s]HGT2.tsv` | Identical open reading frames shared between contigs that have different taxonomic assignments. |

Navigate using the sidebar.
"""
)
