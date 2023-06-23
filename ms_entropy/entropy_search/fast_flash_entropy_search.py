import numpy as np

try:
    from .fast_flash_entropy_search_cpython import  cy_entropy_similarity_identity_search as entropy_similarity_search_identity

except ImportError:
    def entropy_similarity_search_identity(
        product_mz_idx_min,
        product_mz_idx_max,
        intensity,
        entropy_similarity,
        library_peaks_intensity,
        library_spec_idx_array,
        search_spectra_idx_min,
        search_spectra_idx_max,
    ):
        """
        The entropy_similarity will be modified in this function.
        search_type is 0: search all spectra.
        search_type is 1: search spectra in the range [search_spectra_idx_min, search_spectra_idx_max).

        Note: the intensity here should be half of the original intensity.
        """
        all_library_spec_idx = library_spec_idx_array[product_mz_idx_min:product_mz_idx_max]
        idx_list = product_mz_idx_min + np.where(np.bitwise_and(all_library_spec_idx >= search_spectra_idx_min,
                                                 all_library_spec_idx < search_spectra_idx_max))[0]

        array_library_spec_idx = library_spec_idx_array[idx_list]
        array_library_peak_intensity = library_peaks_intensity[idx_list]

        array_library_ab = intensity + array_library_peak_intensity
        entropy_similarity[array_library_spec_idx] += (
            array_library_ab * np.log2(array_library_ab) - intensity * np.log2(intensity) - array_library_peak_intensity * np.log2(array_library_peak_intensity)
        )
