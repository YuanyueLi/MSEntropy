import numpy as np
import unittest
import tempfile
import ms_entropy as me


class TestSpectralEntropy(unittest.TestCase):
    def test_clean_spectrum(self):
        peaks = np.array([[41.04, 0.3716], [69.071, 7.917962], [69.071, 100.0], [86.0969, 66.83]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=85, noise_threshold=0.01, min_ms2_difference_in_da=0.05)
        np.testing.assert_allclose(peaks, np.array([[69.071, 1.0]], dtype=np.float32), atol=1e-4)

    def test_spectral_entropy(self):
        peaks = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        peaks = me.clean_spectrum(
            peaks,
            min_mz=0,
            max_mz=1000,
            noise_threshold=0.01,
            min_ms2_difference_in_da=0.05,
            max_peak_num=None,
            normalize_intensity=True,
        )
        entropy = me.calculate_spectral_entropy(peaks)
        self.assertAlmostEqual(entropy, 0.2605222463607788, places=4)


class TestEntropySimilarity(unittest.TestCase):
    def test_unweighted_entropy_similarity(self):
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        similarity = me.calculate_unweighted_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_da=0.05)
        self.assertAlmostEqual(similarity, 0.9826668790176113, places=4)

        similarity = me.calculate_unweighted_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_da=0.015)
        self.assertAlmostEqual(similarity, 0.9780402849658875, places=4)

    def test_unweighted_entropy_similarity_ppm(self):
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        similarity = me.calculate_unweighted_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_ppm=200)
        self.assertAlmostEqual(similarity, 0.9826668790176113, places=4)

        similarity = me.calculate_unweighted_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_ppm=50)
        self.assertAlmostEqual(similarity, 0.9780402849658875, places=4)

    def test_entropy_similarity(self):
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        similarity = me.calculate_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_da=0.05)
        self.assertAlmostEqual(similarity, 0.8984397722577456, places=4)

        similarity = me.calculate_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_da=0.015)
        self.assertAlmostEqual(similarity, 0.8380874603833756, places=4)

    def test_entropy_similarity_ppm(self):
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        similarity = me.calculate_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_ppm=200)
        self.assertAlmostEqual(similarity, 0.8984397722577456, places=4)

        similarity = me.calculate_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_ppm=50)
        self.assertAlmostEqual(similarity, 0.8380874603833756, places=4)


if __name__ == "__main__":
    unittest.main()
