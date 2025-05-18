import streamlit as st
import matplotlib.pyplot as plt
from rdkit.Chem import Draw

from .core import extract_mols_and_legends, compute_scaffolds

def visualize_molecules(df, title="Top Molecules"):
    """
    Displays grid images of molecules with legends showing score/QED.
    """
    st.header(f"🧪 {title}")
    st.write("Preview of selected compounds with dock score and QED: Quantitative Estimate of Drug-likeness.")

    mols, legends = extract_mols_and_legends(df)
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=4,
        legends=legends,
        subImgSize=(200, 200),
        useSVG=True
    )
    st.image(img, use_container_width=True)

def plot_qed_vs_score(df):
    """
    Scatter plot: Docking score vs QED.
    """
    st.subheader("📊 Docking Score vs. QED")
    fig, ax = plt.subplots()
    ax.scatter(df["DockScore"], df["QED"], alpha=0.6)
    ax.set_xlabel("Docking Score")
    ax.set_ylabel("QED")
    ax.set_title("Docking Score vs. QED")
    st.pyplot(fig)

def show_scaffolds(df):
    """
    Show common Murcko scaffolds.
    """
    st.subheader("🔬 Common Murcko Scaffolds")
    mols, legends = compute_scaffolds(df)
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=3,
        legends=legends,
        subImgSize=(200, 200),
        useSVG=True
    )
    st.image(img, use_container_width=True)

def display_hit_table(df):
    """
    Display final filtered dataframe in a Streamlit table.
    """
    st.dataframe(df)

def convert_to_csv(df):
    """
    Convert dataframe to CSV for download.
    """
    return df.to_csv(index=False).encode("utf-8")
