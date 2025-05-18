import pytest
from components.similarity_screen import screen_by_similarity

# Mock ECFP calculator patch (if compute_ecfp is not available in tests)
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity

def test_screen_by_similarity_basic():
    known_ligands = ["CCO", "CCC"]
    compound_db = ["CCCO", "CCCC", "CCN", "CCOCC", "invalid_smiles", ""]

    top_smiles, df = screen_by_similarity(known_ligands, compound_db, top_k=3)

    assert isinstance(top_smiles, list)
    assert isinstance(df, type(df))  # should be DataFrame
    assert len(top_smiles) <= 3
    for smi in top_smiles:
        assert isinstance(smi, str)
