#define SPEC_TYPE
typedef double float_spec;

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

//' @title Calculate spectral entropy of a spectrum
//' @description Calculate spectral entropy of a spectrum
//'
//' @param peaks A matrix of peaks, with two columns: m/z and intensity.
//' @return A double value of spectral entropy.
//'
//' @examples
//' mz <- c(100.212, 300.321, 535.325)
//' intensity <- c(37.16, 66.83, 999.0)
//' peaks <- matrix(c(mz, intensity), ncol = 2, byrow = FALSE)
//' calculate_spectral_entropy(peaks)
//'
// [[Rcpp::export(name = "calculate_spectral_entropy")]]
double r_calculate_spectral_entropy(const Rcpp::NumericMatrix peaks) {
    Rcpp::NumericVector peaks_vec = convert_matrix_to_vector(peaks);
    int n = peaks_vec.size() / 2;
    const double* peaks_ptr = peaks_vec.begin();
    return calculate_spectral_entropy(peaks_ptr, n);
}

//' @title Clean a spectrum
//' @description Clean a spectrum
//'
//' This function will clean the peaks by the following steps:
//' 1. Remove empty peaks (mz <= 0 or intensity <= 0).
//' 2. Remove peaks with mz >= max_mz or mz < min_mz.
//' 3. Centroid the spectrum by merging peaks within min_ms2_difference_in_da or min_ms2_difference_in_ppm.
//' 4. Remove peaks with intensity < noise_threshold * max_intensity.
//' 5. Keep only the top max_peak_num peaks.
//' 6. Normalize the intensity to sum to 1.
//'
//' Note: The only one of min_ms2_difference_in_da and min_ms2_difference_in_ppm should be positive.
//'
//' @param peaks A matrix of spectral peaks, with two columns: mz and intensity
//' @param min_mz The minimum mz value to keep, set to -1 to disable
//' @param max_mz The maximum mz value to keep, set to -1 to disable
//' @param noise_threshold The noise threshold, set to -1 to disable, all peaks have intensity < noise_threshold * max_intensity will be removed
//' @param min_ms2_difference_in_da The minimum mz difference in Da to merge peaks, set to -1 to disable, any two peaks with mz difference < min_ms2_difference_in_da will be merged
//' @param min_ms2_difference_in_ppm The minimum mz difference in ppm to merge peaks, set to -1 to disable, any two peaks with mz difference < min_ms2_difference_in_ppm will be merged
//' @param max_peak_num The maximum number of peaks to keep, set to -1 to disable
//' @param normalize_intensity Whether to normalize the intensity to sum to 1
//'
//' @return A matrix of spectral peaks, with two columns: mz and intensity
//' @export
//'
//' @examples
//' mz <- c(100.212, 169.071, 169.078, 300.321)
//' intensity <- c(0.3716, 7.917962, 100., 66.83)
//' peaks <- matrix(c(mz, intensity), ncol = 2, byrow = FALSE)
//' clean_spectrum(peaks, min_mz = 0, max_mz = 1000, noise_threshold = 0.01,
//'                min_ms2_difference_in_da = 0.02, min_ms2_difference_in_ppm = -1,
//'                max_peak_num = 100, normalize_intensity = TRUE)
//'
// [[Rcpp::export(name = "clean_spectrum")]]
Rcpp::NumericMatrix r_clean_spectrum(const Rcpp::NumericMatrix peaks,
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
    colnames(peaks_mat) = CharacterVector::create("mz", "intensity");
    return peaks_mat;
}

//' @title Unweighted entropy similarity between two spectra
//' @description Calculate the unweighted entropy similarity between two spectra
//'
//'
//' @param peaks_a A matrix of spectral peaks, with two columns: mz and intensity
//' @param peaks_b A matrix of spectral peaks, with two columns: mz and intensity
//' @param ms2_tolerance_in_da The MS2 tolerance in Da, set to -1 to disable
//' @param ms2_tolerance_in_ppm The MS2 tolerance in ppm, set to -1 to disable
//' @param clean_spectra Whether to clean the spectra before calculating the entropy similarity, see \code{\link{clean_spectrum}}
//' @param min_mz The minimum mz value to keep, set to -1 to disable
//' @param max_mz The maximum mz value to keep, set to -1 to disable
//' @param noise_threshold The noise threshold, set to -1 to disable, all peaks have intensity < noise_threshold * max_intensity will be removed
//' @param max_peak_num The maximum number of peaks to keep, set to -1 to disable
//'
//' @return The unweighted entropy similarity
//'
//' @examples
//' mz_a <- c(169.071, 186.066, 186.0769)
//' intensity_a <- c(7.917962, 1.021589, 100.0)
//' mz_b <- c(120.212, 169.071, 186.066)
//' intensity_b <- c(37.16, 66.83, 999.0)
//' peaks_a <- matrix(c(mz_a, intensity_a), ncol = 2, byrow = FALSE)
//' peaks_b <- matrix(c(mz_b, intensity_b), ncol = 2, byrow = FALSE)
//' calculate_unweighted_entropy_similarity(peaks_a, peaks_b,
//'                                        ms2_tolerance_in_da = 0.02, ms2_tolerance_in_ppm = -1,
//'                                        clean_spectra = TRUE, min_mz = 0, max_mz = 1000,
//'                                        noise_threshold = 0.01,
//'                                        max_peak_num = 100)
//'
// [[Rcpp::export(name = "calculate_unweighted_entropy_similarity")]]
double r_calculate_unweighted_entropy_similarity(const Rcpp::NumericMatrix peaks_a,
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

//' @title Entropy similarity between two spectra
//' @description Calculate the entropy similarity between two spectra
//'
//'
//' @param peaks_a A matrix of spectral peaks, with two columns: mz and intensity
//' @param peaks_b A matrix of spectral peaks, with two columns: mz and intensity
//' @param ms2_tolerance_in_da The MS2 tolerance in Da, set to -1 to disable
//' @param ms2_tolerance_in_ppm The MS2 tolerance in ppm, set to -1 to disable
//' @param clean_spectra Whether to clean the spectra before calculating the entropy similarity, see \code{\link{clean_spectrum}}
//' @param min_mz The minimum mz value to keep, set to -1 to disable
//' @param max_mz The maximum mz value to keep, set to -1 to disable
//' @param noise_threshold The noise threshold, set to -1 to disable, all peaks have intensity < noise_threshold * max_intensity will be removed
//' @param max_peak_num The maximum number of peaks to keep, set to -1 to disable
//'
//' @return The entropy similarity
//'
//' @examples
//' mz_a <- c(169.071, 186.066, 186.0769)
//' intensity_a <- c(7.917962, 1.021589, 100.0)
//' mz_b <- c(120.212, 169.071, 186.066)
//' intensity_b <- c(37.16, 66.83, 999.0)
//' peaks_a <- matrix(c(mz_a, intensity_a), ncol = 2, byrow = FALSE)
//' peaks_b <- matrix(c(mz_b, intensity_b), ncol = 2, byrow = FALSE)
//' calculate_entropy_similarity(peaks_a, peaks_b,
//'                              ms2_tolerance_in_da = 0.02, ms2_tolerance_in_ppm = -1,
//'                              clean_spectra = TRUE, min_mz = 0, max_mz = 1000,
//'                              noise_threshold = 0.01,
//'                              max_peak_num = 100)
//'
// [[Rcpp::export(name = "calculate_entropy_similarity")]]
double r_calculate_entropy_similarity(const Rcpp::NumericMatrix peaks_a,
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
