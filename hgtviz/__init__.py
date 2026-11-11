"""hgt_viz: browse kairos-lr HGT calls contig by contig."""

RANKS = ["d", "p", "c", "o", "f", "g", "s"]
RANK_NAMES = {
    "d": "Domain", "p": "Phylum", "c": "Class", "o": "Order",
    "f": "Family", "g": "Genus", "s": "Species",
}

# Display names for the two kairos-lr HGT criteria (pipeline files keep HGT1 / HGT2).
# Change them here; every page reads these.
HGT1_NAME = "Foreign region"        # region taxonomy differs from its own host contig
HGT2_NAME = "Cross-lineage share"   # same region on contigs of different lineages
HGT1_LABEL = f"{HGT1_NAME} (HGT1)"
HGT2_LABEL = f"{HGT2_NAME} (HGT2)"
