import streamlit as st
from .core import filter_by_quantum_props

def run_quantum_filter_ui(smiles_list):
    st.header("⚛️ Quantum Property Filtering")
    st.write("Filtering based on molecular weight and logP...")

    filtered = filter_by_quantum_props(smiles_list)

    st.success(f"{len(filtered):,} molecules passed quantum-inspired filters.")
    return filtered
