from .core import extract_mols_and_legends, compute_scaffolds

try:
    from .ui import (
        visualize_molecules,
        plot_qed_vs_score,
        show_scaffolds,
        display_hit_table,
        convert_to_csv
    )
    __all__ = [
        "extract_mols_and_legends",
        "compute_scaffolds",
        "visualize_molecules",
        "plot_qed_vs_score",
        "show_scaffolds",
        "display_hit_table",
        "convert_to_csv"
    ]
except ImportError:
    __all__ = [
        "extract_mols_and_legends",
        "compute_scaffolds"
    ]
