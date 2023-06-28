import numpy as np
cimport numpy as np

ctypedef np.float32_t float32
ctypedef np.int64_t int_64
ctypedef np.int8_t int_8
ctypedef np.uint32_t uint_32
ctypedef np.int32_t int_32

ctypedef np.float32_t float_spec
ctypedef int bool

cdef extern from "CleanSpectrum.h":
    int clean_spectrum(float_spec* peaks, int peaks_length,
                        float min_mz, float max_mz,
                        float noise_threshold,
                        float min_ms2_difference_in_da, float min_ms2_difference_in_ppm,
                        int max_peak_num,
                        int normalize_intensity);

cdef extern from "SpectralEntropy.h":
    float calculate_unweighted_entropy_similarity(
        float_spec *peaks_a, int peaks_a_len,
        float_spec *peaks_b, int peaks_b_len,
        float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
        bool clean_spectra,
        float min_mz, float max_mz,
        float noise_threshold,
        int max_peak_num);
    float calculate_entropy_similarity(
        float_spec *peaks_a, int peaks_a_len,
        float_spec *peaks_b, int peaks_b_len,
        float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
        bool clean_spectra,
        float min_mz, float max_mz,
        float noise_threshold,
        int max_peak_num);

cpdef np.ndarray[float32, ndim=2] cy_clean_spectrum(
    peaks,
    min_mz = -1,
    max_mz = -1,
    float noise_threshold = 0.01,
    float min_ms2_difference_in_da = 0.05,
    float min_ms2_difference_in_ppm = -1,
    max_peak_num = -1,
    int normalize_intensity = True,
):
    """
    Clean, centroid, and normalize a spectrum with the following steps:

        1. Remove empty peaks (m/z <= 0 or intensity <= 0).
        2. Remove peaks with m/z >= max_mz or m/z < min_mz.
        3. Centroid the spectrum by merging peaks within min_ms2_difference_in_da.
        4. Remove peaks with intensity < noise_threshold * max_intensity.
        5. Keep only the top max_peak_num peaks.
        6. Normalize the intensity to sum to 1.


    Parameters
    ----------
    peaks : np.ndarray in shape (n_peaks, 2), dtype=np.float32 or list[list[float, float]]
        A 2D array of shape (n_peaks, 2) where the first column is m/z and the second column is intensity.

    min_mz : float, optional
        The minimum m/z to keep. Defaults to None, which will skip removing peaks with m/z < min_mz.

    max_mz : float, optional
        The maximum m/z to keep. Defaults to None, which will skip removing peaks with m/z >= max_mz.

    noise_threshold : float, optional
        The minimum intensity to keep. Defaults to 0.01, which will remove peaks with intensity < 0.01 * max_intensity.

    min_ms2_difference_in_da : float, optional
        The minimum m/z difference between two peaks in the resulting spectrum. Defaults to 0.05, which will merge peaks within 0.05 Da. If a negative value is given, the min_ms2_difference_in_ppm will be used instead.

    min_ms2_difference_in_ppm : float, optional
        The minimum m/z difference between two peaks in the resulting spectrum. Defaults to -1, which will use the min_ms2_difference_in_da instead. If a negative value is given, the min_ms2_difference_in_da will be used instead.
        ** Note either min_ms2_difference_in_da or min_ms2_difference_in_ppm must be positive. If both are positive, min_ms2_difference_in_ppm will be used. **

    max_peak_num : int, optional
        The maximum number of peaks to keep. Defaults to None, which will keep all peaks.

    normalize_intensity : bool, optional
        Whether to normalize the intensity to sum to 1. Defaults to True. If False, the intensity will be kept as is.

        _

    Returns
    -------
    np.ndarray in shape (n_peaks, 2), dtype=np.float32
        The cleaned spectrum will be guaranteed to be sorted by m/z in ascending order.

    """
    if len(peaks) == 0:
        return np.zeros((0, 2), dtype=np.float32)
    if min_ms2_difference_in_da < 0 and min_ms2_difference_in_ppm < 0:
        raise ValueError("Either min_ms2_difference_in_da or min_ms2_difference_in_ppm must be positive.")

    cdef np.ndarray[float32, ndim=2] clean_peaks = np.array(peaks, dtype=np.float32, copy=True, order="C")
    min_mz = -1 if min_mz is None else min_mz
    max_mz = -1 if max_mz is None else max_mz
    max_peak_num = -1 if max_peak_num is None else max_peak_num

    cdef int clean_peaks_len = clean_spectrum(
        <float32*>clean_peaks.data,
        clean_peaks.shape[0],
        min_mz,
        max_mz,
        noise_threshold,
        min_ms2_difference_in_da,
        min_ms2_difference_in_ppm,
        max_peak_num,
        normalize_intensity,
    )
    clean_peaks = clean_peaks[:clean_peaks_len]
    return clean_peaks


cpdef float cy_calculate_entropy_similarity(
    peaks_a, peaks_b,
    float ms2_tolerance_in_da = 0.02,
    float ms2_tolerance_in_ppm = -1,
    bool clean_spectra = True,
    float min_mz = -1,
    float max_mz = -1,
    float noise_threshold = 0.01,
    int max_peak_num = -1
):
    """Calculate the entropy similarity between two spectra.

    First, the spectra are cleaned by the `clean_spectrum()` function.

    Then, the entropy based intensity weights are applied to the peaks.

    Finally, the entropy similarity is calculated by the `calculate_unweighted_entropy_similarity()` function.

    The formula for entropy similarity is as follows:

    .. math::
        Similarity = \\frac{1}{2} \\begin{cases}
        0 & \\text{ if } mz_{A,i} \\neq mz_{B,j} \\\\
        \\sum_{i,j} {f(I_{A,i}+I_{B,j}) - f(I_{A,i}) - f(I_{B,j})} & \\text{ if } mz_{A,i} = mz_{B,j}
        \\end{cases}

    .. math::
        \\text{ where } f(x) = x \\log_2(x) \\text{ and } \\sum_{i} I_{A,i} = \\sum_{j} I_{B,j} = 1

    Parameters
    ----------
    peaks_a : np.ndarray in shape (n_peaks, 2), np.float32 or list[list[float, float]]
        The first spectrum to calculate entropy similarity for. The first column is m/z, and the second column is intensity.

    peaks_b : np.ndarray in shape (n_peaks, 2), np.float32 or list[list[float, float]]
        The second spectrum to calculate entropy similarity for. The first column is m/z, and the second column is intensity.

    ms2_tolerance_in_da : float, optional
        The MS2 tolerance in Da. Defaults to 0.02. If this is set to a negative value, ms2_tolerance_in_ppm will be used instead.

    ms2_tolerance_in_ppm : float, optional
        The MS2 tolerance in ppm. Defaults to -1. If this is set to a negative value, ms2_tolerance_in_da will be used instead.

        **Note:** Either `ms2_tolerance_in_da` or `ms2_tolerance_in_ppm` must be positive. If both `ms2_tolerance_in_da` and `ms2_tolerance_in_ppm` are positive, `ms2_tolerance_in_ppm` will be used.

    clean_spectra : bool, optional
        Whether to clean the spectra before calculating entropy similarity. Defaults to True. **Only set this to False if the spectra have been preprocessed by the `clean_spectrum()` function!** Otherwise, the results will be incorrect. If the spectra are already cleaned, set this to False to save time.

    **kwargs : optional
        The arguments and keyword arguments to pass to function ``clean_spectrum()``.

        _

    Returns
    -------
    float
        The entropy similarity between the two spectra.
    """
    cdef np.ndarray[float32, ndim=2] peaks_a_np = np.array(peaks_a, dtype=np.float32, copy=True, order="C")
    cdef np.ndarray[float32, ndim=2] peaks_b_np = np.array(peaks_b, dtype=np.float32, copy=True, order="C")

    cdef float result = calculate_entropy_similarity(
        <float32*> peaks_a_np.data, len(peaks_a_np),
        <float32*> peaks_b_np.data, len(peaks_b_np),
        ms2_tolerance_in_da,
        ms2_tolerance_in_ppm,
        clean_spectra,
        min_mz,
        max_mz,
        noise_threshold,
        max_peak_num
    )
    return result


cpdef float cy_calculate_unweighted_entropy_similarity(
    peaks_a, peaks_b,
    float ms2_tolerance_in_da = 0.02,
    float ms2_tolerance_in_ppm = -1,
    bool clean_spectra = True,
    float min_mz = -1,
    float max_mz = -1,
    float noise_threshold = 0.01,
    int max_peak_num = -1
):
    """Calculate the unweighted entropy similarity between two spectra.

    The formula for unweighted entropy similarity is as follows:

    .. math::
        Similarity = \\frac{1}{2} \\begin{cases}
        0 & \\text{ if } mz_{A,i} \\neq mz_{B,j} \\\\
        \\sum_{i,j} {f(I_{A,i}+I_{B,j}) - f(I_{A,i}) - f(I_{B,j})} & \\text{ if } mz_{A,i} = mz_{B,j}
        \\end{cases}

    .. math::
        \\text{ where } f(x) = x \\log_2(x) \\text{ and } \\sum_{i} I_{A,i} = \\sum_{j} I_{B,j} = 1


    Parameters
    ----------
    peaks_a : np.ndarray in shape (n_peaks, 2), np.float32 or list[list[float, float]]
        The first spectrum to calculate unweighted entropy similarity for. The first column is m/z, and the second column is intensity.

    peaks_b : np.ndarray in shape (n_peaks, 2), np.float32 or list[list[float, float]]
        The second spectrum to calculate unweighted entropy similarity for. The first column is m/z, and the second column is intensity.

    ms2_tolerance_in_da : float, optional
        The MS2 tolerance in Da. Defaults to 0.02. If this is set to a negative value, ms2_tolerance_in_ppm will be used instead.

    ms2_tolerance_in_ppm : float, optional
        The MS2 tolerance in ppm. Defaults to -1. If this is set to a negative value, ms2_tolerance_in_da will be used instead.

        **Note:** Either `ms2_tolerance_in_da` or `ms2_tolerance_in_ppm` must be positive. If both `ms2_tolerance_in_da` and `ms2_tolerance_in_ppm` are positive, `ms2_tolerance_in_ppm` will be used.

    clean_spectra : bool, optional
        Whether to clean the spectra before calculating unweighted entropy similarity. Defaults to True. Only set this to False if the spectra have been preprocessed by the clean_spectrum() function! Otherwise, the results will be incorrect. If the spectra are already cleaned, set this to False to save time. If the spectra are in the list format, always set this to True or an error will be raised.

    **kwargs : optional
        The arguments and keyword arguments to pass to function ``clean_spectrum()``.

        _

    Returns
    -------
    float
        The unweighted entropy similarity between the two spectra.
    """
    cdef np.ndarray[float32, ndim=2] peaks_a_np
    cdef np.ndarray[float32, ndim=2] peaks_b_np
    if clean_spectra:
        peaks_a_np = np.array(peaks_a, dtype=np.float32, copy=True, order="C")
        peaks_b_np = np.array(peaks_b, dtype=np.float32, copy=True, order="C")
    else:
        peaks_a_np = np.array(peaks_a, dtype=np.float32, copy=False, order="C")
        peaks_b_np = np.array(peaks_b, dtype=np.float32, copy=False, order="C")
        
    cdef float result = calculate_unweighted_entropy_similarity(
        <float32*> peaks_a_np.data, len(peaks_a_np),
        <float32*> peaks_b_np.data, len(peaks_b_np),
        ms2_tolerance_in_da,
        ms2_tolerance_in_ppm,
        clean_spectra,
        min_mz,
        max_mz,
        noise_threshold,
        max_peak_num
    )
    return result
