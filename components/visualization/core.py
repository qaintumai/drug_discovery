from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def extract_mols_and_legends(df, top_n=12):
    mols = [Chem.MolFromSmiles(smi) for smi in df["SMILES"].head(top_n)]
    legends = [f"{row.DockScore:.2f}, QED: {row.QED:.2f}" for _, row in df.head(top_n).iterrows()]
    return mols, legends

def compute_scaffolds(df, top_n=6):
    scaffolds = {}
    for smi in df["SMILES"]:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            smi_scf = Chem.MolToSmiles(scaffold)
            scaffolds[smi_scf] = scaffolds.get(smi_scf, 0) + 1

    sorted_scaffolds = sorted(scaffolds.items(), key=lambda x: x[1], reverse=True)[:top_n]
    mols = [Chem.MolFromSmiles(s[0]) for s in sorted_scaffolds]
    legends = [f"{count} hits" for _, count in sorted_scaffolds]
    return mols, legends
