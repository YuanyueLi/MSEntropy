// SpectralEntropy.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#define __DEBUG__CLEAN_SPECTRUM__ false

#define false 0
#define true 1
#define bool int
// typedef int bool;  // or #define bool int
typedef float float_spec;
// static_assert(sizeof(float_spec) == 4);

void print_spectrum(char* info, float_spec (*spectrum_2d)[2], int spectrum_len);
void swap(float_spec* a, float_spec* b);
void swap_int(int* a, int* b);
void sort_spectrum_by_mz(float_spec (*spectrum_2d)[2], int spectrum_len);
void sort_spectrum_by_mz_and_zero_intensity(float_spec (*spectrum_2d)[2], int* spectrum_len);
void calculate_spectrum_argsort(float_spec (*spectrum_2d)[2], int spectrum_len, int* spectrum_argsort);
bool need_centroid(float_spec (*spectrum_2d)[2], int spectrum_len, float ms2_da);
// Centroid the spectrum, the content in the spectrum will be modified.
int centroid_spectrum(float_spec (*spectrum_2d)[2], int* spectrum_length, float_spec ms2_da, int* spectrum_argsort);

// Clean the spectrum.
// The spectrum is a 2D array. spectrum[x][0] is the m/z, spectrum[x][1] is the intensity.
// The cleaned spectrum will stay the same location as input spectrum.
// The spectrum and spectrum_len will be overwrite by the cleaned spectrum.
// If min_mz, max_mz, noise_threshold, max_peak_num is setted to -1, the corresponding function will be disabled.
int clean_spectrum(float_spec* spectrum, int* spectrum_length,
                   float min_mz, float max_mz,
                   float noise_threshold,
                   int max_peak_num,
                   bool normalize_intensity,
                   float ms2_da);