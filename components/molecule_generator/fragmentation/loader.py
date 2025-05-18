import os
from typing import List
from rdkit import Chem

def load_fragment_library(file_path: str) -> List[str]:
    """
    Load a curated fragment library from a text file of SMILES.

    Args:
        file_path (str): Path to the file containing fragment SMILES, one per line

    Returns:
        List[str]: List of valid fragment SMILES strings
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Fragment library file not found: {file_path}")

    fragments = []
    with open(file_path, 'r') as f:
        for line in f:
            smi = line.strip()
            if smi:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    fragments.append(Chem.MolToSmiles(mol))  # canonicalize
                else:
                    print(f"Warning: Invalid SMILES skipped: {smi}")

    return fragments
