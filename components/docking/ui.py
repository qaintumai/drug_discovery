# components/docking/ui.py

import streamlit as st
from .core import mock_docking

def run_docking_with_ui(smiles_list):
    """
    Run docking with visual feedback using Streamlit UI.
    """
    st.header("🧪 Docking Simulation")
    st.write("Scoring molecules via mock docking...")

    results = mock_docking(smiles_list)

    st.success(f"Docking completed for {len(results):,} molecules.")
    return results