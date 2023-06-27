import numpy as np
from typing import Union
from . import tools


def calculate_entropy_similarity(
    peaks_a,
    peaks_b,
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
        Similarity = \\frac{1}{2} \\sum_{i,j} \\begin{cases}
        0 & \\text{ if } mz_{A,i} \\neq mz_{B,j} \\\\
        f(I_{A,i}+I_{B,j}) - f(I_{A,i}) - f(I_{B,j}) & \\text{ if } mz_{A,i} = mz_{B,j}
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
    if clean_spectra:
        kwargs.update(
            {
                "min_ms2_difference_in_da": max(2 * ms2_tolerance_in_da, kwargs.get("min_ms2_difference_in_da", -1)),
                "min_ms2_difference_in_ppm": max(2 * ms2_tolerance_in_ppm, kwargs.get("min_ms2_difference_in_ppm", -1)),
            }
        )
        peaks_a = tools.clean_spectrum(peaks_a, **kwargs)
        peaks_b = tools.clean_spectrum(peaks_b, **kwargs)
    else:
        peaks_a = np.asarray(peaks_a, dtype=np.float32, order="C").reshape(-1, 2)
        peaks_b = np.asarray(peaks_b, dtype=np.float32, order="C").reshape(-1, 2)

    # Apply the weights to the peaks.
    peaks_a = apply_weight_to_intensity(peaks_a)
    peaks_b = apply_weight_to_intensity(peaks_b)

    return calculate_unweighted_entropy_similarity(
        peaks_a, peaks_b, ms2_tolerance_in_da=ms2_tolerance_in_da, ms2_tolerance_in_ppm=ms2_tolerance_in_ppm, clean_spectra=False
    )


def calculate_unweighted_entropy_similarity(
    peaks_a,
    peaks_b,
    ms2_tolerance_in_da: float = 0.02,
    ms2_tolerance_in_ppm: float = -1,
    clean_spectra: bool = True,
    **kwargs
):
    """Calculate the unweighted entropy similarity between two spectra.

    The formula for unweighted entropy similarity is as follows:

    .. math::
        Similarity = \\frac{1}{2} \\sum_{i,j} \\begin{cases}
        0 & \\text{ if } mz_{A,i} \\neq mz_{B,j} \\\\
        f(I_{A,i}+I_{B,j}) - f(I_{A,i}) - f(I_{B,j}) & \\text{ if } mz_{A,i} = mz_{B,j}
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
        kwargs.update(
            {
                "min_ms2_difference_in_da": max(2 * ms2_tolerance_in_da, kwargs.get("min_ms2_difference_in_da", -1)),
                "min_ms2_difference_in_ppm": max(2 * ms2_tolerance_in_ppm, kwargs.get("min_ms2_difference_in_ppm", -1)),
            }
        )
        peaks_a = tools.clean_spectrum(peaks_a, **kwargs)
        peaks_b = tools.clean_spectrum(peaks_b, **kwargs)
    else:
        peaks_a = np.asarray(peaks_a, dtype=np.float32, order="C").reshape(-1, 2)
        peaks_b = np.asarray(peaks_b, dtype=np.float32, order="C").reshape(-1, 2)

    if peaks_a.shape[0] == 0 or peaks_b.shape[0] == 0:
        return 0.0

    # Calculate the entropy similarity of the two spectra.
    a: int = 0
    b: int = 0
    peak_a_intensity: float = 0.0
    peak_b_intensity: float = 0.0
    peak_ab_intensity: float = 0.0
    entropy_similarity: float = 0.0

    max_allowed_mass_difference: float = ms2_tolerance_in_da

    while a < peaks_a.shape[0] and b < peaks_b.shape[0]:
        mass_difference: float = peaks_a[a, 0] - peaks_b[b, 0]
        if ms2_tolerance_in_ppm > 0:
            max_allowed_mass_difference = peaks_a[a, 0] * ms2_tolerance_in_ppm * 1e-6
        if mass_difference < -max_allowed_mass_difference:
            # This peak only exists in peaks_a.
            a += 1
        elif mass_difference > max_allowed_mass_difference:
            # This peak only exists in peaks_b.
            b += 1
        else:
            # This peak exists in both peaks_a and peaks_b.
            peak_a_intensity = peaks_a[a, 1]
            peak_b_intensity = peaks_b[b, 1]
            peak_ab_intensity = peak_a_intensity + peak_b_intensity
            entropy_similarity += (
                peak_ab_intensity * np.log2(peak_ab_intensity) - peak_a_intensity * np.log2(peak_a_intensity) - peak_b_intensity * np.log2(peak_b_intensity)
            )
            a += 1
            b += 1

    entropy_similarity /= 2
    return entropy_similarity


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


def calculate_spectral_entropy(peaks, clean_spectrum=True, **kwargs) -> float:
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
        peaks = tools.clean_spectrum(peaks, **kwargs)
    else:
        peaks = np.asarray(peaks, dtype=np.float32, order="C").reshape((-1, 2))

    # Calculate the spectral entropy.
    if peaks.shape[0] == 0:
        return 0.0
    else:
        spectral_entropy = -np.sum(peaks[:, 1] * np.log(peaks[:, 1]))
        return spectral_entropy
