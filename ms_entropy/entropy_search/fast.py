import numpy as np


def entropy_similarity_search_fast(
    product_mz_idx_min,
    product_mz_idx_max,
    intensity,
    mixed_spectra_entropy,
    library_peaks_intensity,
    library_spec_idx_array,
    search_type,
    search_spectra_idx_min,
    search_spectra_idx_max,
    search_array,
):
    """
    The mixed_spectra_entropy will be modified in this function.
    search_type is 0: search all spectra.
    search_type is 1: search spectra in the range [search_spectra_idx_min, search_spectra_idx_max).
    search_type is 2: search spectra in the array search_array with entry equals 1, the length of search_array should be equal to the self.total_spectra_num

    Note: the intensity here should be half of the original intensity.
    """
    intensity_xlog2x: np.float32 = intensity * np.log2(intensity)

    for idx in range(product_mz_idx_min, product_mz_idx_max):
        library_spec_idx = library_spec_idx_array[idx]
        if (
            (search_type == 0)
            or (search_type == 1 and search_spectra_idx_min <= library_spec_idx and library_spec_idx < search_spectra_idx_max)
            or (search_type == 2 and search_array[library_spec_idx])
        ):
            # Match this peak
            library_peak_intensity = library_peaks_intensity[idx]
            intensity_ab = intensity + library_peak_intensity

            mixed_spectra_entropy[library_spec_idx] += (
                intensity_ab * np.log2(intensity_ab) - intensity_xlog2x - library_peak_intensity * np.log2(library_peak_intensity)
            )
