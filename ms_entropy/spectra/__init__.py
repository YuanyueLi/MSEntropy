from .tools import clean_spectrum

try:
    from .entropy_numba import calculate_spectral_entropy, calculate_entropy_similarity, calculate_unweighted_entropy_similarity, apply_weight_to_intensity
except ImportError:
    from .entropy import calculate_spectral_entropy, calculate_entropy_similarity, calculate_unweighted_entropy_similarity, apply_weight_to_intensity
