import os
import glob
from pathlib import Path

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


def detect_assembly(data_dir: str) -> str | None:
    matches = glob.glob(os.path.join(data_dir, "*.allHGT1.tsv"))
    if matches:
        return os.path.basename(matches[0]).replace(".allHGT1.tsv", "")
    return None


DEFAULT_DATA_DIR = os.environ.get("HGT_DATA_DIR", str(Path.cwd()))
DEFAULT_ASSEMBLY = os.environ.get("HGT_ASSEMBLY") or detect_assembly(DEFAULT_DATA_DIR) or ""
