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
        Tuple[List[str], pd.DataFrame]: Top-k SMILES and a DataFrame with similarity scores and base ligand.
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

        # Find the max similarity and which known ligand gave it
        sims = [DataStructs.TanimotoSimilarity(fp, qfp) for qfp in query_fps]
        if not sims:
            continue
        max_sim = max(sims)
        max_idx = sims.index(max_sim)
        base_ligand = known_ligands[max_idx]

        results.append((smi, max_sim, base_ligand))

    results.sort(key=lambda x: x[1], reverse=True)
    top_hits = results[:top_k]

    # Return list of SMILES and a dataframe with base ligand info
    top_smiles = [smi for smi, _, _ in top_hits]
    df = pd.DataFrame(top_hits, columns=["SMILES", "Similarity", "Base Ligand"])
    return top_smiles, df
