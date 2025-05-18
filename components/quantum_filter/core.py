from rdkit import Chem
from rdkit.Chem import Descriptors

def filter_by_quantum_props(smiles_list):
    """
    Filter molecules based on simple quantum-inspired descriptors:
    molecular weight and logP.

    Args:
        smiles_list (List[str]): SMILES strings of molecules

    Returns:
        List[str]: SMILES strings that passed the filters
    """
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
