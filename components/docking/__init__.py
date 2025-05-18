# components/docking/__init__.py

# Always available (no streamlit required)
from .core import mock_docking

# Optional: Only expose UI parts if streamlit is available
try:
    from .ui import run_docking_with_ui
    __all__ = ['mock_docking', 'run_docking_with_ui']
except ImportError:
    __all__ = ['mock_docking']