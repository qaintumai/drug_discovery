# tests/__init__.py

import pytest
from pathlib import Path
import sys

# Make sure Python finds components/utils from tests
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.append(str(PROJECT_ROOT))

# Optional: Define fixtures available to all test files
@pytest.fixture
def sample_smiles():
    return ["CCO", "CCCN", "c1ccccc1", "Nc1ccnc2ncc(C(=O)NCC3CCCO3)cc12"]

@pytest.fixture
def large_smiles_list():
    return ["CCO"] * 1000

@pytest.fixture
def empty_smiles_list():
    return []