#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "emscripten.h"
#define __DEBUG__CLEAN_SPECTRUM__ 0
#define __DEBUG__ENTROPY_SIMILARTY__ 0

#define false 0
#define true 1
#define bool int
typedef float float_spec;

/////////////////////////// Helper functions ///////////////////////////
// Print spectrum 2d array
void print_spectrum(const char* info, float_spec (*spectrum_2d)[2], int spectrum_len) {
    printf("%s", info);
    int i;
    for (i = 0; i < spectrum_len; i++) {
        printf("%d\t%f\t%f\n", i, spectrum_2d[i][0], spectrum_2d[i][1]);
    }
}

void swap(float_spec* a, float_spec* b) {
    float_spec c = *a;
    *a = *b;
    *b = c;
}

void inline swap_int(int* a, int* b) {
    int c = *a;
    *a = *b;
    *b = c;
}

// Comparator function for qsort
int compare_by_mz(const void* a, const void* b) {
    float_spec* spectrum_a = (float_spec*)a;
    float_spec* spectrum_b = (float_spec*)b;

    // MZ comparison
    if (spectrum_a[0] < spectrum_b[0]) {
        return -1;
    } else if (spectrum_a[0] > spectrum_b[0]) {
        return 1;
    }

    return 0;
}

void sort_spectrum_by_mz(float_spec (*spectrum_2d)[2], int spectrum_len) {
    // Sort the array using qsort
    qsort(spectrum_2d, spectrum_len, sizeof(float_spec[2]), compare_by_mz);
}

// Comparator function for qsort
int compare_by_mz_with_zero_intensity(const void* a, const void* b) {
    float_spec* spectrum_a = (float_spec*)a;
    float_spec* spectrum_b = (float_spec*)b;

    // Intensity check
    if (spectrum_a[1] > 0 && spectrum_b[1] <= 0) {
        return -1;
    } else if (spectrum_a[1] <= 0 && spectrum_b[1] > 0) {
        return 1;
    }

    // MZ comparison if intensities are the same
    if (spectrum_a[0] < spectrum_b[0]) {
        return -1;
    } else if (spectrum_a[0] > spectrum_b[0]) {
        return 1;
    }

    return 0;
}

int sort_spectrum_by_mz_and_zero_intensity(float_spec (*spectrum_2d)[2], int spectrum_len) {
    // Sort the array using qsort
    qsort(spectrum_2d, spectrum_len, sizeof(float_spec[2]), compare_by_mz_with_zero_intensity);

    // Remove the zero intensity
    while (spectrum_len >= 1 && spectrum_2d[spectrum_len - 1][1] <= 0) {
        spectrum_len--;
    }
    return spectrum_len;
}

int partition(float_spec (*spectrum_2d)[2], int* spectrum_argsort, int low, int high) {
    float pivot = spectrum_2d[spectrum_argsort[high]][1];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (spectrum_2d[spectrum_argsort[j]][1] >= pivot) {
            i++;
            swap_int(&spectrum_argsort[i], &spectrum_argsort[j]);
        }
    }
    swap_int(&spectrum_argsort[i + 1], &spectrum_argsort[high]);
    return (i + 1);
}

void quicksort(float_spec (*spectrum_2d)[2], int* spectrum_argsort, int low, int high) {
    if (low < high) {
        int pi = partition(spectrum_2d, spectrum_argsort, low, high);

        quicksort(spectrum_2d, spectrum_argsort, low, pi - 1);
        quicksort(spectrum_2d, spectrum_argsort, pi + 1, high);
    }
}

void inline calculate_spectrum_argsort(float_spec (*spectrum_2d)[2], int spectrum_len, int* spectrum_argsort) {
    for (int i = 0; i < spectrum_len; i++) {
        spectrum_argsort[i] = i;
    }

    quicksort(spectrum_2d, spectrum_argsort, 0, spectrum_len - 1);
}

/////////////////////////// Spectrum processing functions ///////////////////////////
bool inline need_centroid(float_spec (*spectrum_2d)[2], int spectrum_len, float min_ms2_difference_in_da, float min_ms2_difference_in_ppm) {
    for (int i = 0; i < spectrum_len - 1; i++) {
        if (min_ms2_difference_in_ppm > 0) {
            min_ms2_difference_in_da = spectrum_2d[i + 1][0] * min_ms2_difference_in_ppm * 1e-6;
        }
        if (spectrum_2d[i + 1][0] - spectrum_2d[i][0] < min_ms2_difference_in_da) {
            return true;
        }
    }
    return false;
}

// Centroid the spectrum, the content in the spectrum will be modified.
int centroid_spectrum(float_spec (*spectrum_2d)[2], int spectrum_length, float min_ms2_difference_in_da, float min_ms2_difference_in_ppm, int* spectrum_argsort) {
    // Calculate the argsort of the spectrum by intensity.
    calculate_spectrum_argsort(spectrum_2d, spectrum_length, spectrum_argsort);

    // Centroid the spectrum.
    float mz_delta_allowed_left = min_ms2_difference_in_da;
    float mz_delta_allowed_right = min_ms2_difference_in_da;

    for (int i = 0; i < spectrum_length; i++) {
        int idx = spectrum_argsort[i];
        if (min_ms2_difference_in_ppm > 0) {
            mz_delta_allowed_left = spectrum_2d[idx][0] * min_ms2_difference_in_ppm * 1e-6;
            mz_delta_allowed_right = spectrum_2d[idx][0] / (1 - min_ms2_difference_in_ppm * 1e-6);
        }
        if (spectrum_2d[idx][1] > 0) {
            // Find left board for current peak
            int idx_left = idx - 1;
            while (idx_left >= 0 && spectrum_2d[idx][0] - spectrum_2d[idx_left][0] <= mz_delta_allowed_left) {
                idx_left--;
            }

            // Find right board for current peak
            int idx_right = idx + 1;
            while (idx_right < spectrum_length && spectrum_2d[idx_right][0] - spectrum_2d[idx][0] <= mz_delta_allowed_right) {
                idx_right++;
            }

            // Merge the peaks in the board
            float_spec intensity_sum = 0;
            float_spec intensity_weighted_sum = 0;
            for (int i = idx_left + 1; i < idx_right; i++) {
                intensity_sum += spectrum_2d[i][1];
                intensity_weighted_sum += spectrum_2d[i][1] * spectrum_2d[i][0];
                spectrum_2d[i][1] = 0;
            }

            // Write the new peak into the output spectrum
            spectrum_2d[idx][0] = intensity_weighted_sum / intensity_sum;
            spectrum_2d[idx][1] = intensity_sum;
        }
    }
    spectrum_length = sort_spectrum_by_mz_and_zero_intensity(spectrum_2d, spectrum_length);
    return spectrum_length;
}

// Clean the spectrum.
// The spectrum is a 2D array. spectrum[x][0] is the m/z, spectrum[x][1] is the intensity.
// The spectrum will be rewritten.
EMSCRIPTEN_KEEPALIVE int clean_spectrum(float_spec* spectrum, int spectrum_length,
                                        float min_mz, float max_mz,
                                        float noise_threshold,
                                        float min_ms2_difference_in_da, float min_ms2_difference_in_ppm,
                                        int max_peak_num,
                                        bool normalize_intensity) {
    float_spec(*spectrum_2d)[2] = (float_spec(*)[2]) & spectrum[0];
    int* spectrum_argsort = (int*)malloc(spectrum_length * sizeof(int));

    if (__DEBUG__CLEAN_SPECTRUM__) print_spectrum("Input:\n", spectrum_2d, spectrum_length);
    // 1. Remove the peaks by m/z.
    if (min_mz < 0) {
        min_mz = 0;
    }
    float_spec* spectrum_ptr = spectrum;
    float_spec* spectrum_end = spectrum + spectrum_length * 2;
    for (; spectrum_ptr < spectrum_end; spectrum_ptr += 2) {
        if (*spectrum_ptr <= min_mz || (max_mz > 0 && *spectrum_ptr >= max_mz)) {
            spectrum_ptr[1] = 0;
        }
    }
    spectrum_length = sort_spectrum_by_mz_and_zero_intensity(spectrum_2d, spectrum_length);

    if (__DEBUG__CLEAN_SPECTRUM__) print_spectrum("Remove the peaks by m/z:\n", spectrum_2d, spectrum_length);

    // 2. Centroid the spectrum.
    while (need_centroid(spectrum_2d, spectrum_length, min_ms2_difference_in_da, min_ms2_difference_in_ppm)) {
        spectrum_length = centroid_spectrum(spectrum_2d, spectrum_length, min_ms2_difference_in_da, min_ms2_difference_in_ppm, spectrum_argsort);
        if (__DEBUG__CLEAN_SPECTRUM__) print_spectrum("Centroid the spectrum:\n", spectrum_2d, spectrum_length);
    }
    // 3. Remove the peaks with intensity less than the noise_threshold * maximum(intensity).
    if (noise_threshold > 0) {
        float_spec max_intensity = 0;
        for (int i = 0; i < spectrum_length; i++) {
            if (spectrum_2d[i][1] > max_intensity) {
                max_intensity = spectrum_2d[i][1];
            }
        }
        float_spec noise_threshold_intensity = noise_threshold * max_intensity;
        if (__DEBUG__CLEAN_SPECTRUM__) printf("Remove the peaks with intensity less than %f * %f = %f:\n", noise_threshold, max_intensity, noise_threshold_intensity);
        for (int i = 0; i < spectrum_length; i++) {
            if (spectrum_2d[i][1] < noise_threshold_intensity) {
                spectrum_2d[i][1] = 0;
            }
        }
        if (__DEBUG__CLEAN_SPECTRUM__) print_spectrum("Remove the noise:\n", spectrum_2d, spectrum_length);
    }

    // 4. Select top K peaks.
    if (max_peak_num > 0 && max_peak_num < spectrum_length) {
        calculate_spectrum_argsort(spectrum_2d, spectrum_length, spectrum_argsort);
        for (int i = max_peak_num; i < spectrum_length; i++) {
            spectrum_2d[spectrum_argsort[i]][1] = 0;
        }
        if (__DEBUG__CLEAN_SPECTRUM__) print_spectrum("Select top K peaks:\n", spectrum_2d, spectrum_length);
    }
    spectrum_length = sort_spectrum_by_mz_and_zero_intensity(spectrum_2d, spectrum_length);

    // 5. Normalize the intensity to sum to 1.
    if (normalize_intensity) {
        float_spec sum_intensity = 0;
        float_spec* spectrum_ptr = spectrum + 1;
        const float_spec* spectrum_end = spectrum + spectrum_length * 2;
        for (; spectrum_ptr < spectrum_end; spectrum_ptr += 2) {
            sum_intensity += *spectrum_ptr;
        }
        if (sum_intensity > 0) {
            for (spectrum_ptr = spectrum + 1; spectrum_ptr < spectrum_end; spectrum_ptr += 2) {
                *spectrum_ptr /= sum_intensity;
            }
        }
        if (__DEBUG__CLEAN_SPECTRUM__) print_spectrum("Normalize the intensity:\n", spectrum_2d, spectrum_length);
    }

    free(spectrum_argsort);
    return spectrum_length;
}

/////////////////////////// Spectral Entropy functions ///////////////////////////

// Calculate spectral entropy of a peaks. The peaks intensity need to be prenormalized.
EMSCRIPTEN_KEEPALIVE float_spec calculate_spectral_entropy(const float_spec* peaks, int peaks_length) {
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

// Calculate unweighted entropy similarity
EMSCRIPTEN_KEEPALIVE float calculate_unweighted_entropy_similarity(
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
EMSCRIPTEN_KEEPALIVE float calculate_entropy_similarity(
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
