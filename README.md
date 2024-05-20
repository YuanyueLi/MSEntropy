[![DOI](https://zenodo.org/badge/232434019.svg)](https://zenodo.org/badge/latestdoi/232434019)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7972082.svg)](https://doi.org/10.5281/zenodo.7972082)
[![Test MS Entropy package](https://github.com/YuanyueLi/MSEntropy/actions/workflows/run_test.yml/badge.svg?branch=main)](https://github.com/YuanyueLi/MSEntropy/actions/workflows/run_test.yml)
[![Documentation Status](https://readthedocs.org/projects/msentropy/badge/?version=latest)](https://msentropy.readthedocs.io/en/latest/?badge=latest)


If you have any questions, feel free to send me E-mails: mail@yli.one. If you find this package useful, please consider citing the following papers:

> Li, Y., Fiehn, O. **Flash entropy search to query all mass spectral libraries in real time**, _Nat Methods_ **20**, 1475-1478 (2023). [https://doi.org/10.1038/s41592-023-02012-9](https://doi.org/10.1038/s41592-023-02012-9)
> 
> Li, Y., Kind, T., Folz, J. _et al._ **Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification**, _Nat Methods_ **18**, 1524-1531 (2021). [https://doi.org/10.1038/s41592-021-01331-z](https://doi.org/10.1038/s41592-021-01331-z)

# Theoritical Background

`Spectral entropy` is an useful property to measure the complexity of a spectrum. It is inspried by the concept of Shannon entropy in information theory. [(ref)](https://doi.org/10.1038/s41592-021-01331-z)

`Entropy similarity`, which measured spectral similarity based on spectral entropy, has been shown to outperform dot product similarity in compound identification. [(ref)](https://doi.org/10.1038/s41592-021-01331-z)

The calculation of entropy similarity can be accelerated by using the `Flash Entropy Search` algorithm. [(ref)](https://doi.org/10.1038/s41592-023-02012-9)

# How to use this package

This repository contains the source code to calculate spectral entropy and entropy similarity in various programming languages. Also implemented the Flash Entropy Search algorithm in Python.

## For Python users

A detailed tutorial is available here: [https://msentropy.readthedocs.io](https://msentropy.readthedocs.io)

### Installation

```bash
pip install ms_entropy
```

### Usage of Classical entropy functions

```python
import numpy as np
import ms_entropy as me

peaks_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype = np.float32)
peaks_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype = np.float32)

# Calculate entropy similarity.
entropy = me.calculate_spectral_entropy(peaks_query, clean_spectrum = True, min_ms2_difference_in_da = 0.05)
print(f"Spectral entropy is {entropy}.")

# Calculate unweighted entropy similarity.
unweighted_similarity = me.calculate_unweighted_entropy_similarity(peaks_query, peaks_reference, ms2_tolerance_in_da = 0.05)
print(f"Unweighted entropy similarity: {unweighted_similarity}.")

# Calculate entropy similarity.
similarity = me.calculate_entropy_similarity(peaks_query, peaks_reference, ms2_tolerance_in_da = 0.05)
print(f"Entropy similarity: {similarity}.")
```

### Usage of Flash Entropy Search

```python
from ms_entropy import FlashEntropySearch
entropy_search = FlashEntropySearch()
entropy_search.build_index(spectral_library)
entropy_similarity = entropy_search.search(
    precursor_mz=query_spectrum_precursor_mz, peaks=query_spectrum_peaks)
```

## For R users

A document is available here: [https://cran.r-project.org/web/packages/msentropy/msentropy.pdf](https://cran.r-project.org/web/packages/msentropy/msentropy.pdf)

### Installation

```R
install.packages("msentropy")
```

### Usage

```R
library(msentropy)

# Peaks A
mz_a <- c(169.071, 186.066, 186.0769)
intensity_a <- c(7.917962, 1.021589, 100.0)
peaks_a <- matrix(c(mz_a, intensity_a), ncol = 2, byrow = FALSE)

# Peaks B
mz_b <- c(120.212, 169.071, 186.066)
intensity_b <- c(37.16, 66.83, 999.0)
peaks_b <- matrix(c(mz_b, intensity_b), ncol = 2, byrow = FALSE)

# Calculate spectral entropy
spectral_entropy_a <- calculate_spectral_entropy(clean_spectrum(peaks_a, min_ms2_difference_in_da = 0.02))
spectral_entropy_b <- calculate_spectral_entropy(clean_spectrum(peaks_b, min_ms2_difference_in_da = 0.02))

# Calculate entropy similarity
entropy_similarity <- calculate_entropy_similarity(peaks_a, peaks_b, ms2_tolerance_in_da = 0.02)
```

## For C/C++ users

### Usage

```C
#include "SpectralEntropy.h"

// Calculate spectral entropy
{
    int spec_a_len = 3;
    float spec_a[3][2] = {{169.071, 7.917962}, {186.066, 1.021589}, {186.0769, 100.0}};
    
    // The parameters for clean_spectrum function
    int normalize_intensity = 1;
    float ms2_tolerance_in_da = 0.02, ms2_tolerance_in_ppm = -1;
    float min_mz= -1, max_mz = -1;
    float noise_threshold = 0.01;
    int max_peak_num = -1;

    // Alway clean the spectrum before calculating spectral entropy
    spec_a_len = clean_spectrum(*spec_a, spec_a_len, min_mz, max_mz, noise_threshold, max_peak_num, ms2_tolerance_in_da, ms2_tolerance_in_ppm, max_peak_num, normalize_intensity);

    // Calculate spectral entropy
    float spectral_entropy = calculate_spectral_entropy(*spec_a, spec_a_len);

    printf("Spectral Entropy: %f\n", spectral_entropy);
}

// Calculate entropy similarity
{
    int spec_a_len = 3;
    float spec_a[3][2] = {{169.071, 7.917962}, {186.066, 1.021589}, {186.0769, 100.0}};

    int spec_b_len = 3;
    float spec_b[3][2] = {{120.212, 37.16}, {169.071, 66.83}, {186.066, 999.0}};

    // The parameters for calculate_entropy_similarity function.
    int clean_spectra = 1;
    float ms2_tolerance_in_da = 0.02, ms2_tolerance_in_ppm = -1;
    float min_mz= -1, max_mz = -1;
    float noise_threshold = 0.01;
    int max_peak_num = -1;

    // Calculate entropy similarity, the data in spec_a and spec_b will modified.
    float similarity = calculate_entropy_similarity(*spec_a, spec_a_len, *spec_b, spec_b_len, ms2_tolerance_in_da, ms2_tolerance_in_ppm, clean_spectra, min_mz, max_mz, noise_threshold, max_peak_num);
    printf("Entropy Similarity: %f\n", similarity);
}
```

An example is available in folder [languages/c folder](https://github.com/YuanyueLi/MSEntropy/tree/main/languages/c) and [Example.c](https://github.com/YuanyueLi/MSEntropy/blob/main/languages/c/Example.c), [CMakeLists.txt](https://github.com/YuanyueLi/MSEntropy/blob/main/languages/c/CMakeLists.txt)

## For JavaScript users

An example is available in folder [languages/javascript folder](https://github.com/YuanyueLi/MSEntropy/tree/main/languages/javascript) and [example.js](https://github.com/YuanyueLi/MSEntropy/blob/main/languages/javascript/example.js)

Also, refer to [MSViewer repository](https://github.com/YuanyueLi/MSViewer) for a working example of using this package in a web application.
