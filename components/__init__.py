# components/__init__.py

# Docking
from .docking.core import mock_docking
from .docking.ui import run_docking_with_ui

# Input Section
from .input_section import parse_target_input, input_target, fetch_known_ligands

# Molecule Generation
from .molecule_generator.core import generate_molecules
from .molecule_generator.ui import run_generator_ui

# Molecule Generator: Fragmentation
from .molecule_generator.fragmentation import (
    brics_fragment,
    recap_fragment,
    load_fragment_library,
    filter_fragments_by_size,
    save_fragments_to_file
)

# Prioritization
from .prioritization.core import compute_qed_sascore, prioritize
from .prioritization.ui import run_prioritization_ui

# Quantum Filter
from .quantum_filter.core import filter_by_quantum_props
from .quantum_filter.ui import run_quantum_filter_ui

# Similarity Screening
from .similarity_screen.core import screen_by_similarity
from .similarity_screen.ui import run_similarity_screen_ui

# Visualization
from .visualization.core import extract_mols_and_legends, compute_scaffolds
from .visualization.ui import (
    visualize_molecules,
    plot_qed_vs_score,
    show_scaffolds,
    display_hit_table,
    convert_to_csv
)

# Optionally, define __all__ to specify what gets imported with *
__all__ = [
    # Docking
    'mock_docking',
    'run_docking_with_ui',

    # Input
    'parse_target_input',
    'input_target',
    'fetch_known_ligands',

    # Molecule Generation
    'generate_molecules',
    'run_generator_ui',
    'brics_fragment',
    'recap_fragment',
    'load_fragment_library',
    'filter_fragments_by_size',
    'save_fragments_to_file',

    # Prioritization
    'compute_qed_sascore',
    'prioritize',
    'run_prioritization_ui',

    # Quantum Filter
    'filter_by_quantum_props',
    'run_quantum_filter_ui',

    # Similarity Screening
    'screen_by_similarity',
    'run_similarity_screen_ui',

    # Visualization
    'extract_mols_and_legends',
    'compute_scaffolds',
    'visualize_molecules',
    'plot_qed_vs_score',
    'show_scaffolds',
    'display_hit_table',
    'convert_to_csv'
]