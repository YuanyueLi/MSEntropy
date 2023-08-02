#' @title Calculate spectral entropy similarity between two spectra
#'
#' @description
#'
#' `msentropy_similarity` calculates the spectral entropy between two spectra
#' (Li et al. 2021). It is a wrapper function defining defaults for parameters
#' and calling the [calculate_entropy_similarity()] or
#' [calculate_unweighted_entropy_similarity()] functions to perform the
#' calculation.
#'
#' @param peaks_a A two-column numeric matrix with the m/z and intensity values
#'     for peaks of one spectrum.
#'
#' @param peaks_b A two-column numeric matrix with the m/z and intensity values
#'     for peaks of one spectrum.
#' 
#' @param ms2_tolerance_in_da The MS2 tolerance in Da, set to -1 to disable.
#'     Defaults to `ms2_tolerance_in_da = 0.02`.
#' 
#' @param ms2_tolerance_in_ppm The MS2 tolerance in ppm, set to -1 to disable.
#'     Defaults to `ms2_tolerance_in_ppm = -1`.
#' 
#' @param clean_spectra Whether to clean the spectra before calculating the
#'     entropy similarity, see [clean_spectrum()].
#' 
#' @param min_mz The minimum mz value to keep, set to -1 to disable. Defaults to
#'     `min_mz = 0`.
#' 
#' @param max_mz The maximum mz value to keep, set to -1 to disable. Defaults to
#'     `max_mz = 1000`.
#' 
#' @param noise_threshold The noise threshold, set to -1 to disable, all peaks
#'     have intensity < noise_threshold * max_intensity will be removed.
#'     Defaults to `noise_threshold = 0.01`, thus, by default, all peaks with
#'     an intensity less than 1% of the maximum intensity of a spectrum will
#'     be removed.
#' 
#' @param max_peak_num The maximum number of peaks to keep, set to -1 to
#'     disable. Defaults to `max_peak_num = 1000`.
#'
#' @param weighted `logical(1)` whether the weighted or unweighted entropy
#'     similarity should be calculated. Defaults to `weighted = TRUE`, thus
#'     [calculate_entropy_similarity()] is used for the calculation. For
#'     `weighted = FALSE` [calculate_unweighted_entropy_similarity()] is used
#'     instead.
#' 
#' @param ... Optional additional parameters (currently ignored)
#'
#' @return The entropy similarity
#'
#' @references
#'
#' Li, Y., Kind, T., Folz, J. et al. (2021) Spectral entropy outperforms MS/MS
#' dot product similarity for small-molecule compound identification.
#' Nat Methods 18, 1524-1531.
#' \doi{10.1038/s41592-021-01331-z}.
#' 
#' @examples
#' 
#' peaks_a <- cbind(mz = c(169.071, 186.066, 186.0769),
#'     intensity = c(7.917962, 1.021589, 100.0))
#' peaks_b <- cbind(mz = c(120.212, 169.071, 186.066),
#'     intensity <- c(37.16, 66.83, 999.0))
#' msentropy_similarity(peaks_a, peaks_b, ms2_tolerance_in_da = 0.02)
#' @md
msentropy_similarity <- function(peaks_a, peaks_b, ms2_tolerance_in_da = 0.02,
                                 ms2_tolerance_in_ppm = -1,
                                 clean_spectra = TRUE,
                                 min_mz = 0, max_mz = 1000,
                                 noise_threshold = 0.01,
                                 max_peak_num = 100, weighted = TRUE, ...) {
    if (weighted)
        calculate_entropy_similarity(
            peaks_a, peaks_b, ms2_tolerance_in_da = ms2_tolerance_in_da,
            ms2_tolerance_in_ppm = ms2_tolerance_in_ppm,
            clean_spectra = clean_spectra,
            min_mz = min_mz, max_mz = max_mz, noise_threshold = noise_threshold,
            max_peak_num = max_peak_num)
    else
        calculate_unweighted_entropy_similarity(
            peaks_a, peaks_b, ms2_tolerance_in_da = ms2_tolerance_in_da,
            ms2_tolerance_in_ppm = ms2_tolerance_in_ppm,
            clean_spectra = clean_spectra,
            min_mz = min_mz, max_mz = max_mz, noise_threshold = noise_threshold,
            max_peak_num = max_peak_num)
}
