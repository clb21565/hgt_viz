from pathlib import Path

DATA_DIR = str(Path(__file__).parent.parent / "damhusaen_as_hgt")
ASSEMBLY = "damhusaen_as_rp3_assembly"

TAX_LEVELS = ["d", "p", "c", "o", "f", "g", "s"]
TAX_LEVEL_NAMES = {
    "d": "Domain",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species",
}
