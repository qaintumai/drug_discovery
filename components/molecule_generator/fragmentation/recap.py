from rdkit import Chem
from rdkit.Chem import Recap

def recap_fragment(mol: Chem.Mol) -> list[str]:
    """
    Perform RECAP fragmentation on an RDKit molecule.

    Args:
        mol (Chem.Mol): Input molecule

    Returns:
        List[str]: List of SMILES strings for RECAP fragments
    """
    try:
        recap_tree = Recap.RecapDecompose(mol)
        # Get the leaves (terminal fragments)
        leaves = recap_tree.GetLeaves()
        return [Chem.MolToSmiles(leaf) for leaf in leaves]
    except Exception as e:
        print(f"Error in RECAP fragmentation: {e}")
        return []
