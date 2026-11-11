"""What the HGT call types mean and how to read them."""
import streamlit as st

from hgtviz import HGT1_LABEL, HGT1_NAME, HGT2_LABEL, HGT2_NAME

st.title("HGT calls explained")

# Diagrams: SVG drawn with st.image; mid-tone fills and text read on light and dark backgrounds.
HOST_A = "#2a78d6"    # contig of lineage A
HOST_B = "#eb6834"    # contig / DNA of lineage B
NEUTRAL = "#a3a29c"
TEXT = "#8a8984"      # mid grey: readable on light and dark page backgrounds


def _svg(body, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 {h}" width="640" height="{h}" '
            f'font-family="sans-serif" font-size="15">{body}</svg>')


def _contig(x, y, w, color, label):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="14" rx="3" fill="{color}"/>'
            f'<text x="{x}" y="{y - 14}" fill="{TEXT}">{label}</text>')


st.markdown(f"""
kairos-lr first finds **kairos regions**: stretches of a contig whose genes (ORFs) also
occur, at ≥ 99 % identity, on another contig in the catalogue. Neighbouring shared ORFs
closer than `merge_distance` (1 kb by default) are merged into one region. A region is
only a *candidate*: it says "this DNA is on more than one contig", not yet "this DNA
moved between lineages". Two separate tests turn candidates into HGT calls.
""")

c = st.columns(2)
with c[0]:
    st.subheader(HGT1_LABEL)
    st.image(_svg(
        _contig(20, 40, 560, HOST_A, "host contig: classified as lineage A from all its genes")
        + f'<rect x="230" y="34" width="150" height="26" rx="3" fill="{HOST_B}" opacity="0.85"/>'
        + f'<text x="230" y="82" fill="{TEXT}">kairos region: its own genes say lineage B</text>',
        95), width="stretch")
    st.markdown(f"""
The region's genes are classified on their own (mmseqs taxonomy of the extracted region)
and the answer **disagrees with the contig that carries it**. The region looks foreign to
its host.

- Needs one contig only; no partner is required.
- *Detected at phylum* means region and host agree down to domain and first differ at
  phylum. Coarser (domain, phylum) is a stronger signal than genus or species.
- The region classification must pass the `fragments` / `returned` / `agreed` thresholds
  in `get_hgt-v2.py`, so tiny or ambiguous regions are not called.
- Weak spots: short regions classify poorly, and genes common to many taxa (transposases,
  integrases) can be assigned to whichever lineage dominates the database.
""")

with c[1]:
    st.subheader(HGT2_LABEL)
    st.image(_svg(
        _contig(20, 30, 420, HOST_A, "contig 1: lineage A")
        + _contig(180, 110, 440, HOST_B, "contig 2: lineage B")
        + f'<rect x="200" y="24" width="120" height="26" rx="3" fill="{NEUTRAL}"/>'
        + f'<rect x="300" y="104" width="120" height="26" rx="3" fill="{NEUTRAL}"/>'
        + f'<polygon points="200,50 320,50 420,104 300,104" fill="{NEUTRAL}" opacity="0.35"/>'
        + f'<text x="440" y="70" fill="{TEXT}">same ORFs (≥ 99 % id)</text>',
        140), width="stretch")
    st.markdown(f"""
**Two contigs from different lineages share the same region.** Here the test compares the
two *host contigs'* classifications; the shared DNA itself is not classified.

- Needs at least two contigs; one call per pair of regions.
- *Lineages differ at order* means the two host contigs agree down to class and differ at
  order.
- Has no direction: it does not say which contig gave and which received.
- Weak spots: elements that are everywhere (IS elements, common plasmids) link many
  contigs at once; check how many partners a region has and whether the shared part is
  just one or two genes.
""")

st.subheader("How they combine")
st.markdown(f"""
| You see | Meaning |
|---|---|
| {HGT1_NAME} only | the region looks foreign, but no contig of the matching lineage carries it in this catalogue |
| {HGT2_NAME} only | the region is shared across lineages, but its own genes do not clearly point away from the host |
| both | foreign-looking region that is also found on contigs of another lineage: the strongest static evidence |
| **donor–recipient** | a {HGT1_NAME.lower()} whose partner contig belongs to the lineage the region was classified as. That partner is the putative **donor**, the host the putative **recipient** |
| **HGT group** | regions linked through shared ORFs ({HGT2_NAME.lower()} links), joined transitively. Large groups often come from one widespread gene |
| **in situ** | time-series evidence from the e2e `call` stage: the recipient carries the element only after the donor was seen. The only call type that uses timing |

**Gene-map colours (“HGT evidence”).** A gene is coloured by the coarsest rank at which any
evidence touches it: its region's {HGT1_NAME.lower()} call, or a {HGT2_NAME.lower()} segment
that overlaps it. Darker blue = deeper disagreement (domain, phylum); light grey = in a
kairos region without a call; very light grey = no evidence. Use the **Evidence** switch in
the contig browser to show one type at a time.

Pipeline file names keep the original terms (`allHGT1.tsv`, `allHGT2.tsv`); the display
names are set in `hgtviz/__init__.py`.
""")
