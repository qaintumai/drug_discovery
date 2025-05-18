# molecule_generator/core.py

import numpy as np
from rdkit import Chem
from typing import List, Set
import streamlit as st
from pathlib import Path

@st.cache_data(show_spinner=False)
def load_fragments_from_file(filepath: str = "curated_fragments.smi") -> List[str]:
    """
    Load SMILES fragments from a .smi file.

    Args:
        filepath (str): Path to the SMILES file.

    Returns:
        List[str]: List of fragment SMILES strings.
    """
    path = Path(filepath)
    if not path.exists():
        st.error(f"Fragment file not found: {filepath}")
        return []

    with open(path, "r") as f:
        fragments = [line.strip().split()[0] for line in f if line.strip()]
    return fragments

@st.cache_data(show_spinner=True)
def generate_molecules(n: int = 10000, seed: int = 42) -> List[str]:
    """
    Generate unique molecules by combining curated fragments.

    Args:
        n (int): Number of molecules to generate.
        seed (int): Random seed for reproducibility.

    Returns:
        List[str]: Unique, valid SMILES strings.
    """
    fragments = load_fragments_from_file()
    if not fragments:
        return []

    unique_smiles: Set[str] = set()
    np.random.seed(seed)

    max_attempts = 5 * n
    attempts = 0

    while len(unique_smiles) < n and attempts < max_attempts:
        smi = ''.join(np.random.choice(fragments, size=np.random.randint(2, 5)))
        mol = Chem.MolFromSmiles(smi)
        if mol:
            canon_smi = Chem.MolToSmiles(mol)
            unique_smiles.add(canon_smi)
        attempts += 1

    return list(unique_smiles)
