from rdkit import Chem
from rdkit.Chem import QED, Descriptors
import pandas as pd
from typing import List, Tuple

def compute_qed_sascore(smiles: str) -> Tuple[float, float]:
    mol = Chem.MolFromSmiles(smiles)
    qed_score = QED.qed(mol)
    sa_score = Descriptors.NumRotatableBonds(mol)  # Simple proxy
    return qed_score, sa_score

def prioritize(mol_scores: List[Tuple[str, float]]) -> pd.DataFrame:
    ranked = []
    for smi, dock_score in mol_scores:
        qed_score, sa_score = compute_qed_sascore(smi)
        ranked.append((smi, dock_score, qed_score, sa_score))

    ranked.sort(key=lambda x: (-x[2], x[1], x[3]))  # High QED, better dock, low SA
    top_hits = ranked[:500]

    return pd.DataFrame(top_hits, columns=["SMILES", "DockScore", "QED", "SA"])
