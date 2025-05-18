import pytest
from components.quantum_filter import filter_by_quantum_props

def test_filter_by_quantum_props():
    smiles_list = [
        "CCO",            # Ethanol (too small)
        "CC(=O)O",        # Acetic acid (too small)
        "CCN(CC)CC",      # Triethylamine (valid)
        "CCCCCCCCCCCCCCCCCCCC",  # Too big
        "c1ccccc1",       # Benzene (borderline)
    ]
    filtered = filter_by_quantum_props(smiles_list)

    assert isinstance(filtered, list)
    for smi in filtered:
        assert smi in smiles_list
