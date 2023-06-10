#include <stdio.h>

#include "CleanSpectrum.h"
#include "SpectralEntropy.h"

int main() {
    {
        float_spec spectrum_2d[7][2] = {{41.04, 0.3716}, {0., 0.3716}, {69.070, 7.917962}, {69.070, -7.917962}, {69.071, 100.}, {86.0969, 66.83}, {86.01, 10}};
        int spectrum_len = 7;
        int spectrum_output_len = 0;
        print_spectrum("Origin spectrum:\n", spectrum_2d, spectrum_len);
        spectrum_len = clean_spectrum(*spectrum_2d, spectrum_len, 0, -1, 0.01, 0.05, -1, 5, true);
        print_spectrum("Final spectrum:\n", spectrum_2d, spectrum_len);
    }
    {
        float_spec spec_query[3][2] = {{69.071, 7.917962}, {86.066, 1.021589}, {86.0969, 100.0}};
        float_spec spec_reference[3][2] = {{41.04, 37.16}, {69.07, 66.83}, {86.1, 999.0}};
        int spec_query_len = 3, spec_reference_len = 3;
        float similarity = calculate_entropy_similarity(*spec_query, spec_query_len, *spec_reference, spec_reference_len, 0.02, -1, true, -1, -1, 0.01, -1);
        printf("Entropy Similarity: %f\n", similarity);
    }
    {
        float_spec spec_query[3][2] = {{69.071, 7.917962}, {86.066, 1.021589}, {86.0969, 100.0}};
        float_spec spec_reference[3][2] = {{41.04, 37.16}, {69.07, 66.83}, {86.1, 999.0}};
        int spec_query_len = 3, spec_reference_len = 3;
        float similarity = calculate_unweighted_entropy_similarity(*spec_query, spec_query_len, *spec_reference, spec_reference_len, 0.02, -1, true, -1, -1, 0.01, -1);
        printf("Unweighted entropy similarity: %f\n", similarity);
    }
    return 0;
}
