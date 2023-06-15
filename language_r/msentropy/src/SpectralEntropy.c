#include "SpectralEntropy.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CleanSpectrum.h"

// Calculate unweighted entropy similarity
float calculate_unweighted_entropy_similarity(
    float_spec* peaks_a, int peaks_a_len,
    float_spec* peaks_b, int peaks_b_len,
    float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
    bool clean_spectra,
    float min_mz, float max_mz,
    float noise_threshold,
    int max_peak_num) {
    float_spec(*spec_a_2d)[2] = (float_spec(*)[2])peaks_a;
    float_spec(*spec_b_2d)[2] = (float_spec(*)[2])peaks_b;

    if (__DEBUG__ENTROPY_SIMILARTY__) {
        print_spectrum("spec_query:\n", spec_a_2d, peaks_a_len);
        print_spectrum("spec_reference:\n", spec_b_2d, peaks_b_len);
    }

    if (clean_spectra) {
        peaks_a_len = clean_spectrum(peaks_a, peaks_a_len, min_mz, max_mz, noise_threshold, 2 * ms2_tolerance_in_da, 2 * ms2_tolerance_in_ppm, max_peak_num, true);
        peaks_b_len = clean_spectrum(peaks_b, peaks_b_len, min_mz, max_mz, noise_threshold, 2 * ms2_tolerance_in_da, 2 * ms2_tolerance_in_ppm, max_peak_num, true);
    }

    if (peaks_a_len == 0 || peaks_b_len == 0) {
        return 0.0;
    }

    int a = 0, b = 0;
    float_spec peak_a_intensity, peak_b_intensity, peak_ab_intensity;
    float_spec similarity = 0;

    while (a < peaks_a_len && b < peaks_b_len) {
        float mass_delta_da = spec_a_2d[a][0] - spec_b_2d[b][0];
        if (ms2_tolerance_in_ppm > 0) {
            ms2_tolerance_in_da = ms2_tolerance_in_ppm * spec_a_2d[a][0] * 1e-6;
        }
        if (mass_delta_da < -ms2_tolerance_in_da) {
            // Peak only existed in spec a.
            a++;
        } else if (mass_delta_da > ms2_tolerance_in_da) {
            // Peak only existed in spec b.
            b++;
        } else {
            // Peak existed in both spec a and spec b.
            peak_a_intensity = spec_a_2d[a][1];
            peak_b_intensity = spec_b_2d[b][1];
            peak_ab_intensity = peak_a_intensity + peak_b_intensity;
            similarity += peak_ab_intensity * log2f(peak_ab_intensity) - peak_a_intensity * log2f(peak_a_intensity) - peak_b_intensity * log2f(peak_b_intensity);
            a++;
            b++;
        }
    }
    return similarity / 2;
}

// Calculate entropy similarity
float calculate_entropy_similarity(
    float_spec* peaks_a, int peaks_a_len,
    float_spec* peaks_b, int peaks_b_len,
    float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
    bool clean_spectra,
    float min_mz, float max_mz,
    float noise_threshold,
    int max_peak_num) {
    if (__DEBUG__ENTROPY_SIMILARTY__) {
        print_spectrum("spec_query:\n", (float_spec(*)[2])peaks_a, peaks_a_len);
        print_spectrum("spec_reference:\n", (float_spec(*)[2])peaks_b, peaks_b_len);
    }

    if (clean_spectra) {
        peaks_a_len = clean_spectrum(peaks_a, peaks_a_len, min_mz, max_mz, noise_threshold, 2 * ms2_tolerance_in_da, 2 * ms2_tolerance_in_ppm, max_peak_num, true);
        peaks_b_len = clean_spectrum(peaks_b, peaks_b_len, min_mz, max_mz, noise_threshold, 2 * ms2_tolerance_in_da, 2 * ms2_tolerance_in_ppm, max_peak_num, true);
    }

    if (peaks_a_len == 0 || peaks_b_len == 0) {
        return 0.0;
    }
    apply_weight_to_intensity(peaks_a, peaks_a_len);
    apply_weight_to_intensity(peaks_b, peaks_b_len);

    return calculate_unweighted_entropy_similarity(
        peaks_a, peaks_a_len, peaks_b, peaks_b_len,
        ms2_tolerance_in_da, ms2_tolerance_in_ppm, false, min_mz, max_mz, noise_threshold, max_peak_num);
}

// Calculate spectral entropy of a peaks. The peaks intensity need to be prenormalized.
float_spec calculate_spectral_entropy(const float_spec* peaks, int peaks_length) {
    const float_spec* peak_ptr = &peaks[1];
    const float_spec* peak_end_ptr = peak_ptr + peaks_length * 2;

    float_spec intensity_sum = 0;
    for (; peak_ptr < peak_end_ptr; peak_ptr += 2) {
        if (*peak_ptr > 0) {
            intensity_sum += *peak_ptr;
        }
    }
    if (intensity_sum == 0) {
        return 0;
    } else {
        float_spec entropy = 0;
        for (peak_ptr = &peaks[1]; peak_ptr < peak_end_ptr; peak_ptr += 2) {
            if (*peak_ptr > 0) {
                float_spec intensity = (*peak_ptr) / intensity_sum;
                entropy -= intensity * logf(intensity);
            }
        }
        return entropy;
    }
}

// Apply weight to a peaks by spectral entropy.
// The peaks intensity need to be prenormalized.
// The peaks data will be modified.
void apply_weight_to_intensity(float_spec* peaks, int peaks_length) {
    float_spec entropy = calculate_spectral_entropy(peaks, peaks_length);
    if (entropy < 3) {
        const float_spec weight = 0.25 + 0.25 * entropy;
        float_spec* peak_ptr = &peaks[1];
        const float_spec* peak_end_ptr = peak_ptr + peaks_length * 2;

        // Calculate the sum of intensity.
        float_spec intensity_sum = 0;
        for (; peak_ptr < peak_end_ptr; peak_ptr += 2) {
            *peak_ptr = powf(*peak_ptr, weight);
            intensity_sum += *peak_ptr;
        }

        // Normalize the intensity.
        if (intensity_sum > 0) {
            float_spec reciprocal_intensity_sum = 1.0 / intensity_sum;
            peak_ptr = &peaks[1];
            for (; peak_ptr < peak_end_ptr; peak_ptr += 2) {
                *peak_ptr *= reciprocal_intensity_sum;
            }
        }
    }
}
