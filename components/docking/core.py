# components/docking/core.py

import numpy as np
from typing import List, Tuple

def mock_docking(smiles_list: List[str]) -> List[Tuple[str, float]]:
    """
    Simulates docking by assigning random scores to input molecules.

    Args:
        smiles_list (List[str]): List of SMILES strings representing molecules

    Returns:
        List[Tuple[str, float]]: List of (SMILES, score) pairs
    """
    np.random.seed(42)  # Ensure reproducibility
    results = [(smi, np.random.uniform(-10, -5)) for smi in smiles_list]
    return results