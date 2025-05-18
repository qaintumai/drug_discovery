import streamlit as st
from .core import screen_by_similarity

def run_similarity_screen_ui(known_ligands, compound_db, top_k=1000):
    st.header("🔍 Ligand-Based Similarity Screening")

    with st.spinner("Computing ECFP fingerprints and similarity scores..."):
        top_smiles, df = screen_by_similarity(known_ligands, compound_db, top_k=top_k)

    st.success(f"Top {top_k:,} similar compounds selected.")
    st.dataframe(df.head(10))
    return top_smiles, df
