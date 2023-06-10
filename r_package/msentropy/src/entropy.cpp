#include <Rcpp.h>
using namespace Rcpp;
#include "CleanSpectrum.c"
#include "SpectralEntropy.c"

Rcpp::NumericVector convert_matrix_to_vector(const Rcpp::NumericMatrix peaks) {
    Rcpp::NumericVector peaks_vec = Rcpp::NumericVector(peaks.size());
    // Fill the vector.
    float_spec* peaks_vec_ptr = peaks_vec.begin();
    for (int i = 0; i < peaks.nrow(); i++) {
        *peaks_vec_ptr = peaks(i, 0);
        *(peaks_vec_ptr + 1) = peaks(i, 1);
        peaks_vec_ptr += 2;
    }
    return peaks_vec;
}

Rcpp::NumericMatrix convert_vector_to_matrix(const Rcpp::NumericVector peaks, int nrow) {
    Rcpp::NumericMatrix peaks_mat = Rcpp::NumericMatrix(nrow, 2);
    // Fill the matrix.
    const float_spec* peaks_ptr = peaks.begin();
    for (int i = 0; i < peaks_mat.nrow(); i++) {
        peaks_mat(i, 0) = *peaks_ptr;
        peaks_mat(i, 1) = *(peaks_ptr + 1);
        peaks_ptr += 2;
    }
    return peaks_mat;
}

// [[Rcpp::export]]
float_spec spectralEntropy(const Rcpp::NumericMatrix peaks) {
    Rcpp::NumericVector peaks_vec = convert_matrix_to_vector(peaks);
    int n = peaks_vec.size() / 2;
    const double* peaks_ptr = peaks_vec.begin();
    return calculate_spectral_entropy(peaks_ptr, n);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cleanSpectrum(const Rcpp::NumericMatrix peaks,
                                  float min_mz, float max_mz,
                                  float noise_threshold,
                                  float min_ms2_difference_in_da, float min_ms2_difference_in_ppm,
                                  int max_peak_num,
                                  bool normalize_intensity) {
    Rcpp::NumericVector peaks_vec = convert_matrix_to_vector(peaks);
    int peaks_length = peaks_vec.size() / 2;
    double* peaks_ptr = peaks_vec.begin();
    peaks_length = clean_spectrum(peaks_ptr, peaks_length,
                                  min_mz, max_mz,
                                  noise_threshold,
                                  min_ms2_difference_in_da, min_ms2_difference_in_ppm,
                                  max_peak_num,
                                  normalize_intensity);
    Rcpp::NumericMatrix peaks_mat = convert_vector_to_matrix(peaks_vec, peaks_length);
    return peaks_mat;
}

// [[Rcpp::export]]
double unweightedEntropySimilarity(const Rcpp::NumericMatrix peaks_a,
                                   const Rcpp::NumericMatrix peaks_b,
                                   float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
                                   bool clean_spectra,
                                   float min_mz, float max_mz,
                                   float noise_threshold,
                                   int max_peak_num) {
    Rcpp::NumericVector peaks_a_vec = convert_matrix_to_vector(peaks_a);
    int peaks_a_len = peaks_a_vec.size() / 2;
    double* peaks_a_ptr = peaks_a_vec.begin();
    Rcpp::NumericVector peaks_b_vec = convert_matrix_to_vector(peaks_b);
    int peaks_b_len = peaks_b_vec.size() / 2;
    double* peaks_b_ptr = peaks_b_vec.begin();
    return calculate_unweighted_entropy_similarity(
        peaks_a_ptr, peaks_a_len,
        peaks_b_ptr, peaks_b_len,
        ms2_tolerance_in_da, ms2_tolerance_in_ppm,
        clean_spectra,
        min_mz, max_mz,
        noise_threshold,
        max_peak_num);
}

// [[Rcpp::export]]
double entropySimilarity(const Rcpp::NumericMatrix peaks_a,
                         const Rcpp::NumericMatrix peaks_b,
                         float ms2_tolerance_in_da, float ms2_tolerance_in_ppm,
                         bool clean_spectra,
                         float min_mz, float max_mz,
                         float noise_threshold,
                         int max_peak_num) {
    Rcpp::NumericVector peaks_a_vec = convert_matrix_to_vector(peaks_a);
    int peaks_a_len = peaks_a_vec.size() / 2;
    double* peaks_a_ptr = peaks_a_vec.begin();
    Rcpp::NumericVector peaks_b_vec = convert_matrix_to_vector(peaks_b);
    int peaks_b_len = peaks_b_vec.size() / 2;
    double* peaks_b_ptr = peaks_b_vec.begin();
    return calculate_entropy_similarity(
        peaks_a_ptr, peaks_a_len,
        peaks_b_ptr, peaks_b_len,
        ms2_tolerance_in_da, ms2_tolerance_in_ppm,
        clean_spectra,
        min_mz, max_mz,
        noise_threshold,
        max_peak_num);
}