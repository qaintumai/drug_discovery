# tests/test_input_section.py

from components.input_section import parse_target_input, fetch_known_ligands

def test_parse_target_input_valid_query():
    """Test that valid target queries return expected results"""
    result = parse_target_input("EGFR")
    assert result is not None
    assert result["uniprot_id"] == "P00533"
    assert result["gene_symbol"] == "EGFR"
    assert result["chembl_id"] == "CHEMBL203"
    assert result["organism"] == "Homo sapiens"
    assert result["structure_source"] in ["PDB", "AlphaFold"]

def test_parse_target_input_uniprot_id():
    """Test lookup by UniProt ID"""
    result = parse_target_input("P00533")
    assert result is not None
    assert result["uniprot_id"] == "P00533"

def test_parse_target_input_chembl_id():
    """Test lookup by ChEMBL ID"""
    result = parse_target_input("CHEMBL203")
    assert result is not None
    assert result["chembl_id"] == "CHEMBL203"
    assert result["gene_symbol"] == "EGFR"

def test_parse_target_input_invalid():
    """Test invalid input returns None"""
    result = parse_target_input("InvalidTarget123")
    assert result is None

def test_fetch_known_ligands_for_egfr():
    """Test known ligand retrieval for EGFR (ChEMBL203)"""
    ligands = fetch_known_ligands("CHEMBL203")
    assert isinstance(ligands, list)
    assert len(ligands) > 0
    assert all(isinstance(smi, str) for smi in ligands)

def test_fetch_known_ligands_for_brca1():
    """Test known ligand retrieval for BRCA1 (mocked empty)"""
    ligands = fetch_known_ligands("CHEMBL1234567")
    assert isinstance(ligands, list)
    assert len(ligands) == 1
    assert "Cc1cccc(c1Nc2nccc(n2)c3cccnc3)Cl" in ligands

def test_fetch_known_ligands_unknown_target():
    """Test no ligands returned for unknown target"""
    ligands = fetch_known_ligands("CHEMBL999999")
    assert ligands == []