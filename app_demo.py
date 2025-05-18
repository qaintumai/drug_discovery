# app.py
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, Draw, DataStructs
import numpy as np
import pandas as pd
from PIL import Image

st.set_page_config(page_title="Virtual Drug Discovery Pipeline", layout="wide")

st.title("🧬 Drug Discovery Pipeline")
st.markdown("An end-to-end virtual screening and prioritization app.")

# --- Step 0: Input ---
st.header("1️⃣ Input Target Protein and Known Ligands")
target_protein = st.text_input("Enter target protein name or upload structure", "target.pdb")
known_ligands_input = st.text_area("Enter known ligand SMILES (comma-separated)", "CCO, CCCN")
known_ligands = [s.strip() for s in known_ligands_input.split(",") if s.strip()]

# --- Utility Functions ---
def compute_ecfp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

def ligand_similarity_screen(query_ligands, compound_db):
    query_fps = [compute_ecfp(smi) for smi in query_ligands]
    results = []
    for smi in compound_db:
        fp = compute_ecfp(smi)
        max_sim = max(DataStructs.TanimotoSimilarity(fp, qfp) for qfp in query_fps)
        results.append((smi, max_sim))
    results.sort(key=lambda x: x[1], reverse=True)
    return results[:1000]

def generate_molecules(n=10000):
    frags = ["CCN", "c1ccccc1", "C(=O)O", "CC(C)O"]
    mols = []
    for _ in range(n):
        smi = ''.join(np.random.choice(frags, size=np.random.randint(2,5)))
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols.append(Chem.MolToSmiles(mol))
    return list(set(mols))

def mock_docking(smiles_list):
    np.random.seed(0)
    return [(smi, np.random.uniform(-10, -5)) for smi in smiles_list]

def mock_quantum_filter(smiles_list):
    filtered = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        if 150 < mw < 500 and -2 < logp < 5:
            filtered.append(smi)
    return filtered

def compute_qed_sascore(smiles):
    mol = Chem.MolFromSmiles(smiles)
    qed_score = QED.qed(mol)
    sa_score = Descriptors.NumRotatableBonds(mol)
    return qed_score, sa_score

def prioritize(mol_scores):
    ranked = []
    for smi, score in mol_scores:
        qed_score, sa = compute_qed_sascore(smi)
        ranked.append((smi, score, qed_score, sa))
    ranked.sort(key=lambda x: (-x[2], x[1], x[3]))  # High QED, good docking, low SA
    return ranked[:500]

# --- Run Pipeline ---
if st.button("🚀 Run Pipeline"):
    st.subheader("2️⃣ Ligand-Based Screening")
    with st.spinner("Generating compound library..."):
        compound_db = generate_molecules(10000)
    screened = ligand_similarity_screen(known_ligands, compound_db)
    screened_smiles = [smi for smi, _ in screened]
    st.success(f"Retrieved {len(screened_smiles)} top ligand-similar compounds.")

    st.subheader("3️⃣ Docking")
    docked = mock_docking(screened_smiles)
    st.success(f"Docked {len(docked)} compounds. Showing top 5:")
    for smi, score in docked[:5]:
        st.write(f"**{smi}** – Score: {score:.2f}")
        st.image(Draw.MolToImage(Chem.MolFromSmiles(smi)), width=200)

    st.subheader("4️⃣ Quantum Property Filtering")
    quantum_filtered = mock_quantum_filter([smi for smi, _ in docked])
    filtered_docked = [(smi, score) for smi, score in docked if smi in quantum_filtered]
    st.success(f"{len(filtered_docked)} compounds passed property filter.")

    st.subheader("5️⃣ Prioritization")
    final_hits = prioritize(filtered_docked)
    df = pd.DataFrame(final_hits, columns=["SMILES", "DockScore", "QED", "SA"])
    st.dataframe(df.head(10), use_container_width=True)

    # --- Download CSV ---
    st.download_button("📥 Download Final Hits", df.to_csv(index=False), "final_hits.csv", "text/csv")

    st.balloons()
