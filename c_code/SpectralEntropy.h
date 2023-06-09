// SpectralEntropy.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#define __DEBUG__ENTROPY_SIMILARTY__ false
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CleanSpectrum.h"

// Calculate spectral entropy of a spectrum.
// All peaks in the spectrum need to have a positive intensity, or the result will be wrong.
// The spectrum intensity need to be prenormalized to sum equal to 1.
float_spec calculate_spectral_entropy(const float_spec *spectrum, int spectrum_length);

// Apply weight to a spectrum by spectral entropy.
// The spectrum need to be precleaned by clean_spectrum function, or the result will be wrong.
void apply_weight_to_intensity(float_spec *spectrum, int spectrum_length);

// Calculate unweighted entropy similarity, the spectrum_a and spectrum_b will be modified.
// If the spectra_is_preclean is true, the spectrum_a and spectrum_b need to be cleaned before calculating the similarity, or the similarity will be wrong.
float unweighted_entropy_similarity(
    float_spec *spectrum_a, int spectrum_a_len,
    float_spec *spectrum_b, int spectrum_b_len,
    float ms2_da, bool spectra_is_preclean);

float unweighted_entropy_similarity_ori(
    float_spec *spectrum_a, int spectrum_a_len,
    float_spec *spectrum_b, int spectrum_b_len,
    float ms2_da, bool spectra_is_preclean);

// Calculate entropy similarity, the spectrum_a and spectrum_b will be modified.
// If the spectra_is_preclean is true, the spectrum_a and spectrum_b need to be cleaned before calculating the similarity, or the similarity will be wrong.
float entropy_similarity(
    float_spec *spectrum_a, int spectrum_a_len,
    float_spec *spectrum_b, int spectrum_b_len,
    float ms2_da, bool spectra_is_preclean);
