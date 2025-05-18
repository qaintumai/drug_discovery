from .core import screen_by_similarity

try:
    from .ui import run_similarity_screen_ui
    __all__ = ['screen_by_similarity', 'run_similarity_screen_ui']
except ImportError:
    __all__ = ['screen_by_similarity']
