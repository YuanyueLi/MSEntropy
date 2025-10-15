import numpy as np
cimport numpy as np

ctypedef np.float32_t float32
ctypedef np.int64_t int_64
ctypedef np.int8_t int_8
ctypedef np.uint32_t uint_32
from libc.math cimport log2


cpdef void cy_entropy_similarity_identity_search(int_64 product_mz_idx_min, int_64 product_mz_idx_max,
                                                    float32 intensity, float32[:] entropy_similarity,
                                                    const float32[:] library_peaks_intensity, const uint_32[:] library_spec_idx_array,
                                                    int_64 search_spectra_idx_min, int_64 search_spectra_idx_max) noexcept nogil:
    """
    The entropy_similarity will be modified in this function.

    Note: the intensity here should be half of the original intensity.
    """
    cdef uint_32 library_spec_idx
    cdef float32 library_peak_intensity, intensity_ab
    cdef float32 intensity_xlog2x = intensity * log2(intensity)

    for idx in range(product_mz_idx_min, product_mz_idx_max):
        library_spec_idx = library_spec_idx_array[idx]
        if  search_spectra_idx_min <= library_spec_idx and library_spec_idx < search_spectra_idx_max:
            # Match this peak
            library_peak_intensity = library_peaks_intensity[idx]
            intensity_ab = intensity + library_peak_intensity

            entropy_similarity[library_spec_idx] += \
                intensity_ab * log2(intensity_ab) - \
                intensity_xlog2x - \
                library_peak_intensity * log2(library_peak_intensity)
