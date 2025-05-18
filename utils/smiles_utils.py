# utils/smiles_utils.py

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def compute_ecfp(smiles, radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    return None

def is_valid_smiles(smiles):
    """Check if a SMILES string is valid."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

def compute_properties(smiles):
    """Compute molecular properties for a valid SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    return {
        "SMILES": smiles,
        "Molecular Weight": Descriptors.MolWt(mol),
        "logP": Descriptors.MolLogP(mol),
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Ring Count": Descriptors.RingCount(mol),
    }

def passes_lipinski(smiles):
    """Check if a molecule passes Lipinski's Rule of Five."""
    props = compute_properties(smiles)
    if not props:
        return False

    return (
        props["Molecular Weight"] <= 500 and
        props["logP"] <= 5 and
        props["H-Bond Donors"] <= 10 and
        props["H-Bond Acceptors"] <= 10
    )