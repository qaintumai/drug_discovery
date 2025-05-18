import streamlit as st
from .core import prioritize

def run_prioritization_ui(mol_scores):
    st.header("🎯 Prioritization")
    st.write("Ranking by QED (↑), docking score (↓), and synthetic accessibility (↓)...")

    df = prioritize(mol_scores)

    st.success(f"Top {len(df)} candidates selected.")
    st.dataframe(df.head(10))
    return df
