# Always available (no Streamlit required)
from .core import generate_molecules
from .fragmentation import (
    brics_fragment,
    recap_fragment,
    load_fragment_library,
    filter_fragments_by_size,
    save_fragments_to_file
)

# Optional: Only expose UI parts if Streamlit is available
try:
    from .ui import run_generator_ui
    __all__ = [
        'generate_molecules', 'run_generator_ui',
        'brics_fragment', 'recap_fragment',
        'load_fragment_library', 'filter_fragments_by_size', 'save_fragments_to_file'
    ]
except ImportError:
    __all__ = [
        'generate_molecules',
        'brics_fragment', 'recap_fragment',
        'load_fragment_library', 'filter_fragments_by_size', 'save_fragments_to_file'
    ]
