from dataclasses import dataclass, field


@dataclass
class Feature:
    """
    Generic genomic feature for linear diagrams.
    Holds HGT region data now; designed to accept GFF3 entries later
    (CDS, gene, rRNA, etc.) without changing the renderer.
    """
    start: int
    end: int
    feature_type: str   # "hgt1_region" | "hgt2_region" | "gff_cds" | ...
    label: str
    color: str
    strand: int = 0     # +1 forward, -1 reverse, 0 unknown
    tooltip: str = ""
    metadata: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end - self.start
