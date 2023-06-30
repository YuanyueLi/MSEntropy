#include "CleanSpectrum.h"
#include "SpectralEntropy.h"
#include "emscripten.h"

EMSCRIPTEN_KEEPALIVE int test() {
    return 1;
}

EMSCRIPTEN_KEEPALIVE int wasm_clean_spectrum(float_spec* spectrum, int spectrum_length,
                                             float min_mz, float max_mz,
                                             float noise_threshold,
                                             float min_ms2_difference_in_da, float min_ms2_difference_in_ppm,
                                             int max_peak_num,
                                             bool normalize_intensity) {
    return clean_spectrum(spectrum, spectrum_length,
                          min_mz, max_mz,
                          noise_threshold,
                          min_ms2_difference_in_da, min_ms2_difference_in_ppm,
                          max_peak_num,
                          normalize_intensity);
}

EMSCRIPTEN_KEEPALIVE float_spec wasm_calculate_spectral_entropy(const float_spec* peaks, int peaks_length) {
    return calculate_spectral_entropy(peaks, peaks_length);
}

EMSCRIPTEN_KEEPALIVE float wasm_calculate_unweighted_entropy_similarity(
    float_spec* peaks_a, int peaks_a_len,
    float_spec* peaks_b, int peaks_b_len,
    float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
    bool clean_spectra,
    float min_mz, float max_mz,
    float noise_threshold,
    int max_peak_num) {
    return calculate_unweighted_entropy_similarity(
        peaks_a, peaks_a_len,
        peaks_b, peaks_b_len,
        ms2_tolerance_in_da, ms2_tolerance_in_ppm,
        clean_spectra,
        min_mz, max_mz,
        noise_threshold,
        max_peak_num);
}

EMSCRIPTEN_KEEPALIVE float wasm_calculate_entropy_similarity(
    float_spec* peaks_a, int peaks_a_len,
    float_spec* peaks_b, int peaks_b_len,
    float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
    bool clean_spectra,
    float min_mz, float max_mz,
    float noise_threshold,
    int max_peak_num) {
    return calculate_entropy_similarity(
        peaks_a, peaks_a_len,
        peaks_b, peaks_b_len,
        ms2_tolerance_in_da, ms2_tolerance_in_ppm,
        clean_spectra,
        min_mz, max_mz,
        noise_threshold,
        max_peak_num);
}