// SpectralEntropy.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#define __DEBUG__ENTROPY_SIMILARTY__ 0
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CleanSpectrum.h"

/** Calculate spectral entropy of a spectrum.
 *
 * @param spectrum The spectrum to be calculated. A 2D array. spectrum[x][0] is the m/z, spectrum[x][1] is the intensity.
 * All peaks in the spectrum need to have a positive intensity, or the result will be wrong.
 * @param spectrum_length The length of the spectrum.
 *
 * @return The spectral entropy.
 */
float_spec calculate_spectral_entropy(const float_spec *spectrum, int spectrum_length);

/** Calculate unweighted entropy similarity for two spectra.
 *
 * Note: The peaks_a and peaks_b will be modified, if you want to keep the original spectra, please copy them before calling this function.
 * This function will clean the spectra if clean_spectra is true, otherwise, the spectra will be used directly.
 *
 * Only one of min_ms2_difference_in_da and min_ms2_difference_in_ppm should be positive.
 *
 * @param peaks_a The first spectrum.  A 2D array. spectrum[x][0] is the m/z, spectrum[x][1] is the intensity.
 * @param peaks_a_len The length of the first spectrum.
 * @param peaks_b The second spectrum.  A 2D array. spectrum[x][0] is the m/z, spectrum[x][1] is the intensity.
 * @param peaks_b_len The length of the second spectrum.
 * @param ms2_tolerance_in_da The MS2 tolerance in Da, set to -1 if you want to use ppm.
 * @param ms2_tolerance_in_ppm The MS2 tolerance in ppm, set to -1 if you want to use Da.
 * @param clean_spectra Whether to clean the spectra.
 * If set to false, the spectra need to be cleaned before calling this function, or the result will be wrong.
 * If set to false, the following parameters will be ignored.
 * @param min_mz The minimum m/z of the spectra, set to -1 if you want to use the minimum m/z of the spectra.
 * @param max_mz The maximum m/z of the spectra, set to -1 if you want to use the maximum m/z of the spectra.
 * @param noise_threshold The noise threshold, set to -1 if you want to use the noise threshold of the spectra.
 * @param max_peak_num The maximum number of peaks in the spectra, set to -1 if you want to use the maximum number of peaks of the spectra.
 *
 * @return The unweighted entropy similarity of the two spectra.
 */
float calculate_unweighted_entropy_similarity(
    float_spec *peaks_a, int peaks_a_len,
    float_spec *peaks_b, int peaks_b_len,
    float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
    bool clean_spectra,
    float min_mz, float max_mz,
    float noise_threshold,
    int max_peak_num);

/** Calculate entropy similarity for two spectra.
 *
 * Note: The peaks_a and peaks_b will be modified, if you want to keep the original spectra, please copy them before calling this function.
 * This function will clean the spectra if clean_spectra is true, otherwise, the spectra will be used directly.
 *
 * Only one of min_ms2_difference_in_da and min_ms2_difference_in_ppm should be positive.
 *
 * @param peaks_a The first spectrum.  A 2D array. spectrum[x][0] is the m/z, spectrum[x][1] is the intensity.
 * @param peaks_a_len The length of the first spectrum.
 * @param peaks_b The second spectrum.  A 2D array. spectrum[x][0] is the m/z, spectrum[x][1] is the intensity.
 * @param peaks_b_len The length of the second spectrum.
 * @param ms2_tolerance_in_da The MS2 tolerance in Da, set to -1 if you want to use ppm.
 * @param ms2_tolerance_in_ppm The MS2 tolerance in ppm, set to -1 if you want to use Da.
 * @param clean_spectra Whether to clean the spectra before calculating the similarity.
 * If set to false, the spectra need to be cleaned before calling this function, or the result will be wrong.
 * If set to false, the following parameters will be ignored.
 * @param min_mz The minimum m/z of the spectra, -1 means no limit.
 * @param max_mz The maximum m/z of the spectra, -1 means no limit.
 * @param noise_threshold The noise threshold, -1 means no limit.
 * @param max_peak_num The maximum number of peaks in the spectra, -1 means no limit.
 *
 * @return The entropy similarity of the two spectra.
 */
float calculate_entropy_similarity(
    float_spec *peaks_a, int peaks_a_len,
    float_spec *peaks_b, int peaks_b_len,
    float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
    bool clean_spectra,
    float min_mz, float max_mz,
    float noise_threshold,
    int max_peak_num);

// Apply weight to a spectrum by spectral entropy.
// The spectrum need to be precleaned by clean_spectrum function, or the result will be wrong.
void apply_weight_to_intensity(float_spec *spectrum, int spectrum_length);
