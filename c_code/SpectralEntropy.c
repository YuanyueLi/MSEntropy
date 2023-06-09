#include "SpectralEntropy.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CleanSpectrum.h"

// Calculate unweighted entropy similarity
float unweighted_entropy_similarity(
    float_spec* spectrum_a, int spectrum_a_len,
    float_spec* spectrum_b, int spectrum_b_len,
    float ms2_da, bool spectra_is_preclean) {
    if (!spectra_is_preclean) {
        clean_spectrum(spectrum_a, &spectrum_a_len, 0, -1, 0.01, -1, true, 2 * ms2_da);
        clean_spectrum(spectrum_b, &spectrum_b_len, 0, -1, 0.01, -1, true, 2 * ms2_da);
    }

    float_spec(*spec_a_2d)[2] = (float_spec(*)[2])spectrum_a;
    float_spec(*spec_b_2d)[2] = (float_spec(*)[2])spectrum_b;
    if (__DEBUG__ENTROPY_SIMILARTY__) {
        print_spectrum("spec_query:\n", spec_a_2d, spectrum_a_len);
        print_spectrum("spec_reference:\n", spec_b_2d, spectrum_b_len);
    }

    if (spectrum_a_len == 0 || spectrum_b_len == 0) {
        return 0.0;
    }

    int a = 0, b = 0;
    float_spec peak_a_intensity, peak_b_intensity, peak_ab_intensity;
    float_spec similarity = 0;

    while (a < spectrum_a_len && b < spectrum_b_len) {
        float mass_delta_da = spec_a_2d[a][0] - spec_b_2d[b][0];
        if (mass_delta_da < -ms2_da) {
            // Peak only existed in spec a.
            a++;
        } else if (mass_delta_da > ms2_da) {
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
float entropy_similarity(
    float_spec* spectrum_a, int spectrum_a_len,
    float_spec* spectrum_b, int spectrum_b_len,
    float ms2_da, bool spectra_is_preclean) {
    if (!spectra_is_preclean) {
        clean_spectrum(spectrum_a, &spectrum_a_len, 0, -1, 0.01, -1, true, 2 * ms2_da);
        clean_spectrum(spectrum_b, &spectrum_b_len, 0, -1, 0.01, -1, true, 2 * ms2_da);
    }
    if (__DEBUG__ENTROPY_SIMILARTY__) {
        print_spectrum("spec_query:\n", (float_spec(*)[2])spectrum_a, spectrum_a_len);
        print_spectrum("spec_reference:\n", (float_spec(*)[2])spectrum_b, spectrum_b_len);
    }
    if (spectrum_a_len == 0 || spectrum_b_len == 0) {
        return 0.0;
    }
    apply_weight_to_intensity(spectrum_a, spectrum_a_len);
    apply_weight_to_intensity(spectrum_b, spectrum_b_len);

    return unweighted_entropy_similarity(spectrum_a, spectrum_a_len, spectrum_b, spectrum_b_len, ms2_da, true);
}

// Calculate spectral entropy of a spectrum. The spectrum intensity need to be prenormalized.
float_spec calculate_spectral_entropy(const float_spec* spectrum, int spectrum_length) {
    const float_spec* spectrum_ptr = &spectrum[1];
    const float_spec* spectrum_end = spectrum_ptr + spectrum_length * 2;

    float_spec entropy = 0;
    for (; spectrum_ptr < spectrum_end; spectrum_ptr += 2) {
        float_spec intensity = *spectrum_ptr;
        entropy -= intensity * logf(intensity);
    }

    return entropy;
}

// Apply weight to a spectrum by spectral entropy.
// The spectrum intensity need to be prenormalized.
// The spectrum data will be modified.
void apply_weight_to_intensity(float_spec* spectrum, int spectrum_length) {
    float_spec entropy = calculate_spectral_entropy(spectrum, spectrum_length);
    if (entropy < 3) {
        const float_spec weight = 0.25 + 0.25 * entropy;
        float_spec intensity_sum = 0;
        float_spec* spectrum_ptr = &spectrum[1];
        const float_spec* spectrum_end = spectrum_ptr + spectrum_length * 2;
        for (; spectrum_ptr < spectrum_end; spectrum_ptr += 2) {
            *spectrum_ptr = powf(*spectrum_ptr, weight);
            intensity_sum += *spectrum_ptr;
        }
        if (intensity_sum > 0) {
            float_spec reciprocal_intensity_sum = 1.0 / intensity_sum;
            spectrum_ptr = &spectrum[1];
            for (; spectrum_ptr < spectrum_end; spectrum_ptr += 2) {
                *spectrum_ptr *= reciprocal_intensity_sum;
            }
        }
    }
}
