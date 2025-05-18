# utils/chembl.py

import requests

def get_known_ligands(uniprot_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/target?target_components.accession={uniprot_id}&format=json"
    try:
        response = requests.get(url)
        target_data = response.json()
        if target_data["targets"]:
            chembl_id = target_data["targets"][0]["target_chembl_id"]
            activity_url = f"https://www.ebi.ac.uk/chembl/api/data/activity?target_chembl_id={chembl_id}&format=json&limit=1000"
            act_response = requests.get(activity_url)
            activity_data = act_response.json()
            ligands = list({
                act["molecule_chembl_id"]
                for act in activity_data["activities"]
                if act.get("standard_type") in ["IC50", "Ki"]
            })
            # Optional: convert chembl_id to SMILES via ChEMBL (left out for brevity)
            return ligands
    except Exception:
        pass
    return []
