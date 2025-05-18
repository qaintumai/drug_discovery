# app.py
# Streamlit frontend for the drug discovery pipeline

import streamlit as st
from components import (
    input_section,
    similarity_screen,
    molecule_generator,
    docking,
    quantum_filter,
    prioritization,
    visualization,
)

st.set_page_config(page_title="Drug Discovery Pipeline", layout="wide")

# --- App Title ---
st.title("🧬 Drug Discovery Pipeline")
st.markdown("End-to-end virtual screening and prioritization")

# --- Step 0: Target Input ---
st.header("🎯 Step 0: Enter Target Protein")
target_id, known_ligands = input_section.input_target()

if target_id:
    # --- Step 1: Ligand-based Similarity Screening ---
    st.header("🔍 Step 1: Ligand Similarity Screening")

    st.subheader("📌 Known Ligands Retrieved")
    for smi in known_ligands:
        st.code(smi, language="smiles")

    st.info("Generating virtual compound library...")
    compound_db = molecule_generator.generate_molecules(n=10000)

    st.info("Computing Tanimoto similarity...")
    screened_smiles, similarity_df = similarity_screen.screen_by_similarity(known_ligands, compound_db)

    # Optional: Show top 5 with best similarity
    st.subheader("🔍 Top Similar Molecules")
    st.dataframe(similarity_df)  # Will now have columns SMILES, Similarity, Base Ligand

    # --- Step 2: Molecule Generation ---
    st.header("🧪 Step 2: Molecule Generation")
    st.write(f"✅ Generated {len(compound_db)} unique molecules using mock diverse fragments as candidate library.")

    # --- Step 3: Docking ---
    st.header("🧲 Step 3: Docking (Mock Scores)")
    docked = docking.mock_docking(screened_smiles)
    st.success("Docking completed.")

    # --- Step 4: Quantum Property Filtering ---
    st.header("🔬 Step 4: Quantum Property Filtering")
    quantum_pass = quantum_filter.filter_by_quantum_props([smi for smi, _ in docked])
    filtered_docked = [(smi, score) for smi, score in docked if smi in quantum_pass]
    st.success(f"{len(filtered_docked)} compounds passed quantum filtering.")


    # --- Step 5: Prioritization ---
    st.header("📈 Step 5: Prioritization")
    final_hits = prioritization.prioritize(filtered_docked)
    st.success(f"{len(final_hits)} final hits selected for experimental validation.")

    # --- Output Table ---
    st.header("📋 Final Compound Hits")
    visualization.display_hit_table(final_hits)

    st.download_button(
        label="📥 Download Results (CSV)",
        data=visualization.convert_to_csv(final_hits),
        file_name="final_hits.csv",
        mime="text/csv",
    )

    # --- Optional: Molecule Images ---
    visualization.visualize_molecules(final_hits[:20], title="Top 20 Molecules")
    visualization.plot_qed_vs_score(final_hits)
    visualization.show_scaffolds(final_hits)

else:
    st.warning("Please enter a valid protein target or UniProt ID to start the pipeline.")
