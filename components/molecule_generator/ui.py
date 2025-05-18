# molecule_generator/ui.py (or similar)

import streamlit as st
from .core import generate_molecules

def run_generator_ui(default_n: int = 10000):
    """
    Streamlit UI for generating molecules using fragment combinations.

    Args:
        default_n (int): Default number of molecules to generate

    Returns:
        List[str]: Generated SMILES strings
    """
    st.header("🧬 Molecule Generation")
    n = st.number_input(
        "Number of molecules to generate",
        min_value=10, max_value=10000,
        value=default_n, step=10
    )

    if st.button("Generate Molecules"):
        with st.spinner(f"Generating {n:,} candidate molecules..."):
            smiles_list = generate_molecules(n)
        st.success(f"{len(smiles_list):,} molecules generated.")
        return smiles_list

    return []
