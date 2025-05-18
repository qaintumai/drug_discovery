import pandas as pd
from components.visualization import extract_mols_and_legends, compute_scaffolds

def test_extract_mols_and_legends():
    df = pd.DataFrame({
        "SMILES": ["CCO", "CCN", "c1ccccc1"],
        "DockScore": [-7.5, -6.8, -5.9],
        "QED": [0.6, 0.7, 0.8],
    })
    mols, legends = extract_mols_and_legends(df, top_n=2)
    assert len(mols) == 2
    assert len(legends) == 2
    assert all(isinstance(l, str) for l in legends)

def test_compute_scaffolds():
    df = pd.DataFrame({
        "SMILES": ["CCO", "CCO", "CCN", "CCN", "c1ccccc1", "c1ccccc1"],
    })
    mols, legends = compute_scaffolds(df, top_n=2)
    assert len(mols) <= 2
    assert len(legends) == len(mols)
    assert all("hits" in l for l in legends)
