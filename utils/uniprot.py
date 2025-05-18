# utils/uniprot.py

import requests

def resolve_protein_name_or_id(query):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json&size=1"
    try:
        response = requests.get(url)
        data = response.json()
        if data["results"]:
            result = data["results"][0]
            return {
                "uniprot_id": result["primaryAccession"],
                "name": result["proteinDescription"]["recommendedName"]["fullName"]["value"]
            }
    except Exception:
        pass
    return {"uniprot_id": "Unknown", "name": query}
