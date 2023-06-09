import numpy as np
from typing import Union
import ctypes
from pathlib import Path

lib_spectral_entropy = ctypes.cdll.LoadLibrary(Path(__file__).parent / "libSpectralEntropy.so")


def tools_clean_spectrum(
    peaks: Union[list[list[float, float]], np.ndarray],
    min_mz: float = -1,
    max_mz: float = -1,
    noise_threshold: float = 0.01,
    min_ms2_difference_in_da: float = 0.05,
    min_ms2_difference_in_ppm: float = -1,
    max_peak_num: int = -1,
    normalize_intensity: bool = True,
    **kwargs
) -> np.ndarray:
    clean_peaks = np.array(peaks, dtype=np.float32, copy=True, order="C")
    cfunc_clean_spectrum = lib_spectral_entropy.clean_spectrum
    cfunc_clean_spectrum.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, flags="C_CONTIGUOUS"),  # peaks
        ctypes.c_int,  # spectrum_length
        ctypes.c_float,  # min_mz
        ctypes.c_float,  # max_mz
        ctypes.c_float,  # noise_threshold
        ctypes.c_float,  # min_ms2_difference_in_da
        ctypes.c_float,  # min_ms2_difference_in_ppm
        ctypes.c_int,  # max_peak_num
        ctypes.c_bool,  # normalize_intensity
    ]
    cfunc_clean_spectrum.restype = ctypes.c_int

    min_mz = -1 if min_mz is None else min_mz
    max_mz = -1 if max_mz is None else max_mz
    max_peak_num = -1 if max_peak_num is None else max_peak_num

    clean_peaks_len = cfunc_clean_spectrum(
        clean_peaks,
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


def calculate_entropy_similarity(
    peaks_a: Union[list[list[float, float]], np.ndarray],
    peaks_b: Union[list[list[float, float]], np.ndarray],
    ms2_tolerance_in_da: float = 0.02,
    ms2_tolerance_in_ppm: float = -1,
    clean_spectra: bool = True,
    **kwargs
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
    peaks_a = np.array(peaks_a, dtype=np.float32, copy=True, order="C")
    peaks_b = np.array(peaks_b, dtype=np.float32, copy=True, order="C")
    cfunc_calculate_entropy_similarity = lib_spectral_entropy.calculate_entropy_similarity
    cfunc_calculate_entropy_similarity.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, flags="C_CONTIGUOUS"),  # peaks_a
        ctypes.c_int,  # len(peaks_a)
        np.ctypeslib.ndpointer(dtype=np.float32, flags="C_CONTIGUOUS"),  # peaks_b
        ctypes.c_int,  # len(peaks_b)
        ctypes.c_float,  # ms2_tolerance_in_da
        ctypes.c_float,  # ms2_tolerance_in_ppm
        ctypes.c_bool,  # clean_spectra
        ctypes.c_float,  # min_mz
        ctypes.c_float,  # max_mz
        ctypes.c_float,  # noise_threshold
        ctypes.c_int,  # max_peak_num
    ]
    cfunc_calculate_entropy_similarity.restype = ctypes.c_float
    result = cfunc_calculate_entropy_similarity(
        peaks_a,
        len(peaks_a),
        peaks_b,
        len(peaks_b),
        ms2_tolerance_in_da,
        ms2_tolerance_in_ppm,
        clean_spectra,
        kwargs.get("min_mz", -1),
        kwargs.get("max_mz", -1),
        kwargs.get("noise_threshold", -1),
        kwargs.get("max_peak_num", -1),
    )
    return result


def calculate_unweighted_entropy_similarity(
    peaks_a: Union[list[list[float, float]], np.ndarray],
    peaks_b: Union[list[list[float, float]], np.ndarray],
    ms2_tolerance_in_da: float = 0.02,
    ms2_tolerance_in_ppm: float = -1,
    clean_spectra: bool = True,
    **kwargs
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
    if clean_spectra:
        peaks_a = np.array(peaks_a, dtype=np.float32, copy=True, order="C")
        peaks_b = np.array(peaks_b, dtype=np.float32, copy=True, order="C")
    cfunc_calculate_unweighted_entropy_similarity = lib_spectral_entropy.calculate_unweighted_entropy_similarity
    cfunc_calculate_unweighted_entropy_similarity.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, flags="C_CONTIGUOUS"),  # peaks_a
        ctypes.c_int,  # len(peaks_a)
        np.ctypeslib.ndpointer(dtype=np.float32, flags="C_CONTIGUOUS"),  # peaks_b
        ctypes.c_int,  # len(peaks_b)
        ctypes.c_float,  # ms2_tolerance_in_da
        ctypes.c_float,  # ms2_tolerance_in_ppm
        ctypes.c_bool,  # clean_spectra
        ctypes.c_float,  # min_mz
        ctypes.c_float,  # max_mz
        ctypes.c_float,  # noise_threshold
        ctypes.c_int,  # max_peak_num
    ]
    cfunc_calculate_unweighted_entropy_similarity.restype = ctypes.c_float
    result = cfunc_calculate_unweighted_entropy_similarity(
        peaks_a,
        len(peaks_a),
        peaks_b,
        len(peaks_b),
        ms2_tolerance_in_da,
        ms2_tolerance_in_ppm,
        clean_spectra,
        kwargs.get("min_mz", -1),
        kwargs.get("max_mz", -1),
        kwargs.get("noise_threshold", -1),
        kwargs.get("max_peak_num", -1),
    )
    return result


def apply_weight_to_intensity(peaks: np.ndarray) -> np.ndarray:
    """
    Apply a weight to the intensity of a spectrum based on spectral entropy based on the method described in:

    Li, Y., Kind, T., Folz, J. et al. Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. Nat Methods 18, 1524\-1531 (2021). https://doi.org/10.1038/s41592-021-01331-z.

    Parameters
    ----------
    peaks : np.ndarray in shape (n_peaks, 2), np.float32
        The spectrum to apply weight to. The first column is m/z, and the second column is intensity.
        The peaks need to be pre-cleaned.

        _

    Returns
    -------
    np.ndarray in shape (n_peaks, 2), np.float32
        The spectrum with weight applied. The first column is m/z, and the second column is intensity.
        The peaks will be a copy of the input peaks.
    """
    if peaks.shape[0] == 0:
        return np.empty((0, 2), dtype=np.float32)

    # Calculate the spectral entropy.
    entropy = 0.0
    if peaks.shape[0] > 0:
        entropy = -np.sum(peaks[:, 1] * np.log(peaks[:, 1]))

    # Copy the peaks.
    weighted_peaks = peaks.copy()

    # Apply the weight.
    if entropy < 3:
        weight = 0.25 + 0.25 * entropy
        weighted_peaks[:, 1] = np.power(peaks[:, 1], weight)
        intensity_sum = np.sum(weighted_peaks[:, 1])
        weighted_peaks[:, 1] /= intensity_sum

    return weighted_peaks


def calculate_spectral_entropy(peaks: Union[list[list[float, float]], np.ndarray], clean_spectrum=True, **kwargs) -> float:
    """Calculate the spectral entropy of a spectrum.

    Parameters
    ----------
    peaks : np.ndarray in shape (n_peaks, 2), np.float32 or list[list[float, float]]
        The spectrum to calculate spectral entropy for. The first column is m/z, and the second column is intensity.

    clean_spectrum : bool, optional
        Whether to clean the spectrum before calculating spectral entropy. Defaults to True. If the spectrum is already cleaned, set this to False to save time.

    **kwargs : optional
        The arguments and keyword arguments to pass to clean_spectrum().

        _

    Returns
    -------
    float
        The spectral entropy of the spectrum.
    """
    # Clean the spectrum.
    if clean_spectrum:
        peaks = tools_clean_spectrum(peaks, **kwargs)
    else:
        peaks = np.asarray(peaks, dtype=np.float32, order="C").reshape((-1, 2))

    # Calculate the spectral entropy.
    if peaks.shape[0] == 0:
        return 0.0
    else:
        spectral_entropy = -np.sum(peaks[:, 1] * np.log(peaks[:, 1]))
        return spectral_entropy
