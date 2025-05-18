import pytest
import pandas as pd
from components.prioritization import compute_qed_sascore, prioritize

# Example SMILES (simple and valid)
TEST_SMILES = [
    ("CCO", -7.2),     # Ethanol
    ("c1ccccc1", -6.5),  # Benzene
    ("CC(=O)O", -8.0),  # Acetic acid
]

def test_compute_qed_sascore_returns_valid_scores():
    for smiles, _ in TEST_SMILES:
        qed, sa = compute_qed_sascore(smiles)
        assert 0.0 <= qed <= 1.0, f"QED score out of range: {qed}"
        assert isinstance(sa, (int, float)), "SA score is not a number"

def test_prioritize_returns_dataframe():
    df = prioritize(TEST_SMILES)
    assert isinstance(df, pd.DataFrame)
    assert all(col in df.columns for col in ["SMILES", "DockScore", "QED", "SA"])
    assert len(df) <= 500  # should cap at 500 top molecules
    assert df.iloc[0]["QED"] >= df.iloc[-1]["QED"] or len(df) == 1
