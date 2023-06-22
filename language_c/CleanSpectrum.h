// SpectralEntropy.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#define __DEBUG__CLEAN_SPECTRUM__ 0

#include <stdbool.h>
// #define false 0
// #define true 1
// #define bool int
// typedef int bool;  // or #define bool int
#ifdef SPEC_TYPE
#else
typedef float float_spec;
#endif
// static_assert(sizeof(float_spec) == 4);

/**
 * @brief Clean the spectrum.
 *
 * The function will modify the content in the peaks in place and return the length of the cleaned peaks.
 * If you want to keep the original peaks, please copy it before calling this function.
 *
 * This function will clean the peaks by the following steps:
 * 1. Remove empty peaks (m/z <= 0 or intensity <= 0).
 * 2. Remove peaks with m/z >= max_mz or m/z < min_mz.
 * 3. Centroid the spectrum by merging peaks within min_ms2_difference_in_da or min_ms2_difference_in_ppm.
 * 4. Remove peaks with intensity < noise_threshold * max_intensity.
 * 5. Keep only the top max_peak_num peaks.
 * 6. Normalize the intensity to sum to 1.
 *
 * Note: The only one of min_ms2_difference_in_da and min_ms2_difference_in_ppm should be positive.

 * @param peaks The peaks to be cleaned. A 2D array. peaks[x][0] is the m/z, peaks[x][1] is the intensity.
 * @param peaks_length The length of the peaks.
 * @param min_mz The minimum m/z of the peaks. If set to -1, this function will not remove peaks with m/z < min_mz.
 * @param max_mz The maximum m/z of the peaks. If set to -1, this function will not remove peaks with m/z >= max_mz.
 * @param noise_threshold The noise threshold of the peaks. If set to -1, this function will not remove peaks with intensity < noise_threshold * max_intensity.
 * @param min_ms2_difference_in_da The minimum difference in m/z to merge peaks. If set to -1, this function will not centroid the peaks.
 * @param min_ms2_difference_in_ppm The minimum difference in ppm to merge peaks. If set to -1, this function will not centroid the peaks.
 * @param max_peak_num The maximum number of peaks to keep. If set to -1, this function will not remove peaks.
 * @param normalize_intensity Whether to normalize the intensity to sum to 1.
 *
 * @return int The length of the cleaned peaks.
*/
int clean_spectrum(float_spec* peaks, int peaks_length,
                   float min_mz, float max_mz,
                   float noise_threshold,
                   float min_ms2_difference_in_da, float min_ms2_difference_in_ppm,
                   int max_peak_num,
                   bool normalize_intensity);

void print_spectrum(const char* info, float_spec (*spectrum_2d)[2], int spectrum_len);
void swap(float_spec* a, float_spec* b);
void swap_int(int* a, int* b);
void sort_spectrum_by_mz(float_spec (*spectrum_2d)[2], int spectrum_len);
int sort_spectrum_by_mz_and_zero_intensity(float_spec (*spectrum_2d)[2], int spectrum_len);
void calculate_spectrum_argsort(float_spec (*spectrum_2d)[2], int spectrum_len, int* spectrum_argsort);
bool need_centroid(float_spec (*spectrum_2d)[2], int spectrum_len, float min_ms2_difference_in_da, float min_ms2_difference_in_ppm);
// Centroid the spectrum, the content in the spectrum will be modified.
int centroid_spectrum(float_spec (*spectrum_2d)[2], int peaks_length, float min_ms2_difference_in_da, float min_ms2_difference_in_ppm, int* spectrum_argsort);
