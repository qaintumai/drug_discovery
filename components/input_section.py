# components/input_section.py

import requests
import streamlit as st
from typing import Tuple, Optional

def parse_target_input(target_query: str) -> dict:
    """
    Accepts Common Name, UniProt ID, ChEMBL ID, or Gene Symbol
    Returns a dictionary with:
        uniprot_id, chembl_id, gene_symbol, organism, structure_source
    """
    # Mock data for demo purposes
    mock_data = {
        "EGFR": {
            "uniprot_id": "P00533",
            "chembl_id": "CHEMBL203",
            "gene_symbol": "EGFR",
            "organism": "Homo sapiens",
            "structure_source": "PDB"
        },
        "BRCA1": {
            "uniprot_id": "Q12984",
            "chembl_id": "CHEMBL1234567",
            "gene_symbol": "BRCA1",
            "organism": "Homo sapiens",
            "structure_source": "AlphaFold"
        },
        "P00533": {
            "uniprot_id": "P00533",
            "chembl_id": "CHEMBL203",
            "gene_symbol": "EGFR",
            "organism": "Homo sapiens",
            "structure_source": "PDB"
        },
        "CHEMBL203": {
            "uniprot_id": "P00533",
            "chembl_id": "CHEMBL203",
            "gene_symbol": "EGFR",
            "organism": "Homo sapiens",
            "structure_source": "PDB"
        },
        "BRCA1_HUMAN": {
            "uniprot_id": "Q12984",
            "chembl_id": "CHEMBL1234567",
            "gene_symbol": "BRCA1",
            "organism": "Homo sapiens",
            "structure_source": "AlphaFold"
        }
    }

    return mock_data.get(target_query.strip(), None)

def input_target():
    """
    Streamlit UI component to accept target input and return parsed info
    """
    st.text_input("Enter Protein Target (Name, UniProt ID, ChEMBL ID, or Gene Symbol)", key="target_query")
    target_query = st.session_state.target_query.strip()

    if not target_query:
        return None, []

    result = parse_target_input(target_query)

    if result is None:
        st.error("Could not resolve the protein target. Please try another identifier.")
        return None, []

    # Show resolved information
    st.markdown(f"✅ Resolved Target:")
    st.markdown(f"- **UniProt ID:** {result['uniprot_id']}")
    st.markdown(f"- **ChEMBL ID:** {result['chembl_id']}")
    st.markdown(f"- **Gene Symbol:** {result['gene_symbol']}")
    st.markdown(f"- **Organism:** {result['organism']}")
    st.markdown(f"- **Recommended Structure:** {result['structure_source']}")

    # Fetch known ligands (mocked)
    known_ligands = fetch_known_ligands(result["chembl_id"])

    return result, known_ligands

def fetch_known_ligands(chembl_id: str):
    """Mock fetching known ligands from ChEMBL"""
    if chembl_id == "CHEMBL203":
        return ["CCOc1ccccc1C2=Nc3cccnc3NC2", "CN(C)c1ccc(cc1)C2=Cc3ccccc3C2=N"]
    elif chembl_id == "CHEMBL1234567":
        return ["Cc1cccc(c1Nc2nccc(n2)c3cccnc3)Cl"]
    return []