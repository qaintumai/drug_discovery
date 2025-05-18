# tests/test_docking.py

from components.docking import mock_docking

def test_mock_docking_output_format():
    smiles_list = ["CCO", "CCCN", "c1ccccc1"]
    results = mock_docking(smiles_list)

    assert isinstance(results, list)
    for item in results:
        assert isinstance(item, tuple)
        assert len(item) == 2
        assert isinstance(item[0], str)
        assert isinstance(item[1], float)

def test_mock_docking_reproducibility():
    smiles_list = ["CCO", "CCCN", "c1ccccc1"]
    results1 = mock_docking(smiles_list)
    results2 = mock_docking(smiles_list)
    assert results1 == results2

def test_mock_docking_empty_input():
    assert mock_docking([]) == []

def test_mock_docking_large_input():
    results = mock_docking(["CCO"] * 1000)
    assert len(results) == 1000