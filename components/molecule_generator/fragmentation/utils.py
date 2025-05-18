from rdkit import Chem
from typing import List

def filter_fragments_by_size(fragments: List[str], min_atoms: int = 3, max_atoms: int = 30) -> List[str]:
    """
    Filter fragments by number of heavy atoms.

    Args:
        fragments (List[str]): List of fragment SMILES
        min_atoms (int): Minimum number of heavy atoms
        max_atoms (int): Maximum number of heavy atoms

    Returns:
        List[str]: Filtered list of fragment SMILES
    """
    filtered = []
    for smi in fragments:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            num_atoms = mol.GetNumHeavyAtoms()
            if min_atoms <= num_atoms <= max_atoms:
                filtered.append(smi)
    return filtered

def save_fragments_to_file(fragments: List[str], file_path: str):
    """
    Save fragments to a text file, one SMILES per line.

    Args:
        fragments (List[str]): List of fragment SMILES
        file_path (str): Output file path
    """
    with open(file_path, 'w') as f:
        for smi in fragments:
            f.write(smi + '\n')
