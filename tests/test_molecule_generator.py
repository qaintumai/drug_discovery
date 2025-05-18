import pytest
from components.molecule_generator import generate_molecules
from rdkit import Chem

def test_generate_molecule_count():
    num_molecules = 50
    molecules = generate_molecules(n=num_molecules)
    assert isinstance(molecules, list)
    assert len(molecules) == num_molecules

def test_generate_molecules_are_valid_smiles():
    molecules = generate_molecules(n=20)
    for smi in molecules:
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Invalid SMILES: {smi}"
