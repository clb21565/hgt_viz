import os
import pandas as pd
import streamlit as st
from config import DATA_DIR, ASSEMBLY, TAX_LEVELS


def _path(suffix):
    return os.path.join(DATA_DIR, f"{ASSEMBLY}.{suffix}")


@st.cache_data(show_spinner=False)
def load_hgt1():
    df = pd.read_csv(_path("allHGT1.tsv"), sep="\t", low_memory=False)
    # Derive detection level: coarsest level where region taxonomy != contig taxonomy
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
def load_hgt2():
    return pd.read_csv(_path("allHGT2.tsv"), sep="\t", low_memory=False)


@st.cache_data(show_spinner=False)
def load_dr():
    return pd.read_csv(_path("allDR.tsv"), sep="\t", low_memory=False)


@st.cache_data(show_spinner=False)
def load_groups():
    return pd.read_csv(_path("hgt_groups.tsv"), sep="\t", low_memory=False)


@st.cache_data(show_spinner=False)
def load_report():
    df = pd.read_csv(_path("hgt_report.tsv"), sep="\t", low_memory=False)
    # Parse "region_types" string "_internal_:1988,_end_edge_:1430" into separate columns
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
