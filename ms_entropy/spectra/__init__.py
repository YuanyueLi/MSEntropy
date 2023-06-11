from .entropy import calculate_spectral_entropy, apply_weight_to_intensity

try:
    from .entropy_cython import (
        cy_clean_spectrum as clean_spectrum,
        cy_calculate_entropy_similarity as calculate_entropy_similarity,
        cy_calculate_unweighted_entropy_similarity as calculate_unweighted_entropy_similarity,
    )
except ImportError:
    from .tools import clean_spectrum
    from .entropy import calculate_entropy_similarity, calculate_unweighted_entropy_similarity
