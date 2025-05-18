from rdkit import Chem
from rdkit.Chem import BRICS

def brics_fragment(mol: Chem.Mol) -> list[str]:
    """
    Perform BRICS fragmentation on an RDKit molecule.

    Args:
        mol (Chem.Mol): Input molecule

    Returns:
        List[str]: List of SMILES strings for BRICS fragments
    """
    try:
        frags = BRICS.BRICSDecompose(mol)
        return list(frags)
    except Exception as e:
        print(f"Error in BRICS fragmentation: {e}")
        return []
