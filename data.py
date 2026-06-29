import os
import pandas as pd
import streamlit as st
from config import TAX_LEVELS, detect_assembly, DEFAULT_DATA_DIR, DEFAULT_ASSEMBLY


def _path(data_dir: str, assembly: str, suffix: str) -> str:
    return os.path.join(data_dir, f"{assembly}.{suffix}")


@st.cache_data(show_spinner=False)
def load_hgt1(data_dir: str, assembly: str) -> pd.DataFrame:
    df = pd.read_csv(_path(data_dir, assembly, "allHGT1.tsv"), sep="\t", low_memory=False)
    df["detection_level"] = pd.NA
    for lvl in TAX_LEVELS:
        r_col, c_col = f"{lvl}_region", f"{lvl}_contig"
        if r_col in df.columns and c_col in df.columns:
            mask = (
                df["detection_level"].isna()
                & (df[r_col].fillna("") != df[c_col].fillna(""))
                & df[r_col].notna()
                & df[c_col].notna()
            )
            df.loc[mask, "detection_level"] = lvl
    df["detection_level"] = df["detection_level"].fillna("unknown")
    return df


@st.cache_data(show_spinner=False)
def load_hgt2(data_dir: str, assembly: str) -> pd.DataFrame:
    return pd.read_csv(_path(data_dir, assembly, "allHGT2.tsv"), sep="\t", low_memory=False)


@st.cache_data(show_spinner=False)
def load_dr(data_dir: str, assembly: str) -> pd.DataFrame:
    return pd.read_csv(_path(data_dir, assembly, "allDR.tsv"), sep="\t", low_memory=False)


@st.cache_data(show_spinner=False)
def load_groups(data_dir: str, assembly: str) -> pd.DataFrame:
    return pd.read_csv(_path(data_dir, assembly, "hgt_groups.tsv"), sep="\t", low_memory=False)


@st.cache_data(show_spinner=False)
def load_report(data_dir: str, assembly: str) -> pd.DataFrame:
    df = pd.read_csv(_path(data_dir, assembly, "hgt_report.tsv"), sep="\t", low_memory=False)

    def parse_region_types(s):
        result = {"rt_internal": 0, "rt_end_edge": 0, "rt_start_edge": 0}
        if pd.isna(s):
            return result
        for part in str(s).split(","):
            if ":" in part:
                k, v = part.split(":", 1)
                k = k.strip().strip("_")
                col = f"rt_{k.replace(' ', '_')}"
                if col in result:
                    result[col] = int(v.strip())
        return result

    parsed = df["region_types"].apply(parse_region_types).apply(pd.Series)
    return pd.concat([df, parsed], axis=1)


def dataset_sidebar() -> tuple[str, str]:
    """Render dataset controls in the sidebar; return (data_dir, assembly)."""
    with st.sidebar:
        st.header("Dataset")

        if "data_dir" not in st.session_state:
            st.session_state["data_dir"] = DEFAULT_DATA_DIR
        if "assembly" not in st.session_state:
            st.session_state["assembly"] = DEFAULT_ASSEMBLY

        st.text_input("Data directory", key="data_dir")
        data_dir = st.session_state["data_dir"]

        detected = detect_assembly(data_dir)
        # Auto-fill assembly when it's blank (e.g. after user clears it or on first load)
        if detected and not st.session_state.get("assembly"):
            st.session_state["assembly"] = detected

        st.text_input(
            "Assembly prefix",
            key="assembly",
            help="Filename prefix shared by all HGT output files (auto-detected from *.allHGT1.tsv)",
        )
        assembly = st.session_state["assembly"]

        if detected and detected != assembly:
            st.caption(f"Detected in directory: `{detected}`")

    return data_dir, assembly
