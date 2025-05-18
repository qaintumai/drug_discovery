from rdkit import Chem, DataStructs
from utils.smiles_utils import compute_ecfp
import pandas as pd

def screen_by_similarity(known_ligands, compound_db, top_k=1000):
    """
    Perform similarity screening using ECFP fingerprints and Tanimoto similarity.

    Args:
        known_ligands (List[str]): SMILES of known active ligands.
        compound_db (List[str]): SMILES of compounds to screen.
        top_k (int): Number of top results to return.

    Returns:
        Tuple[List[str], pd.DataFrame]: Top-k SMILES and their similarity scores.
    """
    query_fps = [
        compute_ecfp(smi) for smi in known_ligands
        if Chem.MolFromSmiles(smi)
    ]

    results = []
    for smi in compound_db:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        fp = compute_ecfp(smi)
        max_sim = max(
            DataStructs.TanimotoSimilarity(fp, qfp) for qfp in query_fps
        )
        results.append((smi, max_sim))

    results.sort(key=lambda x: x[1], reverse=True)
    top_hits = results[:top_k]

    return [smi for smi, _ in top_hits], pd.DataFrame(top_hits, columns=["SMILES", "Similarity"])
