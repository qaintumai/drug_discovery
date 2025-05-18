# 🧬 Ligand-Based Virtual Screening and Docking App

A streamlined pipeline for ligand-based virtual screening, molecular generation, and docking — designed for drug discovery teams, computational chemists, and bioinformaticians.

Built with **Streamlit**, **RDKit**, and **Python**, this app enables scientists to:

- Input a **target** and retrieve known ligands from ChEMBL
- Screen 10,000+ candidate molecules via **Tanimoto similarity**
- Perform **mock docking** with scored output (placeholder for real docking integration)
- Prioritize molecules based on **QED** (drug-likeness) and **SA** (synthetic accessibility)
- (Optionally) filter candidates using **quantum properties**

---

## 🚀 Features

| Step | Description |
|------|-------------|
| 🔍 1. Ligand Similarity | Retrieve known ligands and compute Tanimoto similarity against a virtual compound library |
| 🧪 2. Molecule Generation | Generate diverse candidate molecules using fragment-based synthesis (BRICS, RECAP) |
| ⚓ 3. Docking | Assign mock docking scores to simulate binding affinity |
| 🔬 4. Prioritization | Score and filter molecules using drug-likeness and synthetic accessibility |
| ⚛️ 5. Quantum Filtering *(optional)* | Apply quantum property constraints (e.g., HOMO-LUMO gap, dipole moment) |

---

## 🛠️ Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
cd YOUR_REPO_NAME
pip install -r requirements.txt
```
## ▶️ Running the App
```bash
streamlit run app.py
```