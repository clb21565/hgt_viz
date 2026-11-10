"""COG functional categories (eggNOG COG_category letters) and the groups used to
colour genes. A gene can carry several letters (e.g. "KT")."""

COG = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}

INFO, CELL, METAB, MOBILE, POOR, NONE = (
    "Information storage & processing", "Cellular processes & signaling", "Metabolism",
    "Mobilome (X)", "Poorly characterized (R, S)", "No COG / no eggNOG hit")
GROUP_OF = {**{c: INFO for c in "JAKLB"}, **{c: CELL for c in "DYVTMNZWUO"},
            **{c: METAB for c in "CGEFHIPQ"}, "X": MOBILE, "R": POOR, "S": POOR}
GROUPS = [MOBILE, INFO, CELL, METAB, POOR, NONE]


def letters(cat) -> list[str]:
    """'KT' -> ['K', 'T']; missing or '-' -> []."""
    if not isinstance(cat, str):
        return []
    return [c for c in cat if c in COG]


def group(cat) -> str:
    """One group per gene: Mobilome if any letter is X (the HGT-relevant call), else
    the group of the first informative letter, else poorly characterized / none."""
    ls = letters(cat)
    if not ls:
        return NONE
    if "X" in ls:
        return MOBILE
    for c in ls:
        if GROUP_OF[c] != POOR:
            return GROUP_OF[c]
    return POOR


def label(letter: str) -> str:
    return f"{letter}: {COG.get(letter, letter)}"
