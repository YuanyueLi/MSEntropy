import numpy as np
import unittest
import ms_entropy as me


class TestSpectralEntropy(unittest.TestCase):
    def test_clean_spectrum(self):
        # Test
        peaks = np.array([[41.04, 0.3716], [69.071, 7.917962], [69.071, 100.0], [86.0969, 66.83]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=85, noise_threshold=0.01, min_ms2_difference_in_da=0.05)
        np.testing.assert_allclose(peaks, np.array([[69.071, 1.0]], dtype=np.float32), atol=1e-4)

        print("*" * 80)
        print("Test start")
        # Test
        peaks = np.array([[41.04, 0.3716], [69.070, 0.917962], [69.071, 100.0], [86.0969, 66.83]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, noise_threshold=0.01, min_ms2_difference_in_ppm=50)
        np.testing.assert_allclose(peaks, np.array([[69.07099, 0.6016047], [86.0969, 0.39839533]], dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([[41.04, 0.3716], [69.070, 0.917962], [69.071, 100.0], [86.0969, 66.83]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=85, noise_threshold=0.01, min_ms2_difference_in_ppm=50)
        np.testing.assert_allclose(peaks, np.array([[69.07099, 1.0]], dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([[41.04, 1.3716], [69.071, 100.0], [86.0969, 66.83]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=85, noise_threshold=0.01, min_ms2_difference_in_ppm=50, normalize_intensity=False)
        np.testing.assert_allclose(peaks, np.array([[41.04, 1.3716], [69.071, 100.0]], dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=85, noise_threshold=0.01, min_ms2_difference_in_da=0.05)
        np.testing.assert_allclose(peaks, np.zeros((0, 2), dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([[41.04, 0.0], [69.071, 0.0], [86.0969, 66.83]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=85, noise_threshold=0.01, min_ms2_difference_in_da=0.05, normalize_intensity=False)
        np.testing.assert_allclose(peaks, np.zeros((0, 2), dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([[86.0969, 66.83]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=85, noise_threshold=0.01, min_ms2_difference_in_da=0.05)
        np.testing.assert_allclose(peaks, np.zeros((0, 2), dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([[41.04, 20.3716], [69.071, 7.917962], [86.0969, 66.83]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=85, min_mz=50, noise_threshold=0.01, min_ms2_difference_in_da=0.05)
        np.testing.assert_allclose(peaks, np.array([[69.071, 1.0]], dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([[41.04, 10], [69.071, 30], [86.0969, 60]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, max_mz=-1, min_mz=-1, noise_threshold=0.01, min_ms2_difference_in_da=0.05)
        np.testing.assert_allclose(peaks, np.array([[41.04, 0.1], [69.071, 0.3], [86.0969, 0.6]], dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([[41.04, 10], [69.071, 40], [86.0969, 60]], dtype=np.float32)
        peaks = me.clean_spectrum(peaks, noise_threshold=0.01, min_ms2_difference_in_da=0.05, max_peak_num=2)
        np.testing.assert_allclose(peaks, np.array([[69.071, 0.4], [86.0969, 0.6]], dtype=np.float32), atol=1e-4)

        # Test
        peaks = np.array([[86.0969, 66.83]], dtype=np.float32)
        self.assertRaises(ValueError, me.clean_spectrum, peaks, noise_threshold=0.01, min_ms2_difference_in_da=-1)

    def test_spectral_entropy(self):
        # Test
        peaks = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        peaks = me.clean_spectrum(
            peaks,
            min_mz=0.0,
            max_mz=1000.0,
            noise_threshold=0.01,
            min_ms2_difference_in_da=0.05,
            max_peak_num=-1,
            normalize_intensity=True,
        )
        entropy = me.calculate_spectral_entropy(peaks, clean_spectrum=False)
        self.assertAlmostEqual(entropy, 0.2605222463607788, places=4)

        # Test
        peaks = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        entropy = me.calculate_spectral_entropy(
            peaks, clean_spectrum=True, min_mz=0, max_mz=1000, noise_threshold=0.01, min_ms2_difference_in_da=0.05, max_peak_num=None, normalize_intensity=True
        )
        self.assertAlmostEqual(entropy, 0.2605222463607788, places=4)

        # Test
        peaks = np.array([[86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        entropy = me.calculate_spectral_entropy(
            peaks, clean_spectrum=True, min_mz=0, max_mz=85, noise_threshold=0.01, min_ms2_difference_in_da=0.05, max_peak_num=None, normalize_intensity=True
        )
        self.assertAlmostEqual(entropy, 0.0, places=4)


class TestEntropySimilarity(unittest.TestCase):
    def test_unweighted_entropy_similarity(self):
        # Test
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        similarity = me.calculate_unweighted_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_da=0.05)
        self.assertAlmostEqual(similarity, 0.9826668790176113, places=4)

        # Test
        similarity = me.calculate_unweighted_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_da=0.015)
        self.assertAlmostEqual(similarity, 0.9780402849658875, places=4)

        # Test
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        similarity = me.calculate_unweighted_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_da=0.05, max_mz=50)
        self.assertAlmostEqual(similarity, 0.0, places=4)

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

        # Test
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        similarity = me.calculate_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_da=0.05, max_mz=50)
        self.assertAlmostEqual(similarity, 0.0, places=4)

        # Test
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        spec_query = me.clean_spectrum(
            spec_query,
            min_mz=0,
            max_mz=85,
            noise_threshold=0.01,
        )
        spec_reference = me.clean_spectrum(
            spec_reference,
            min_mz=50,
            max_mz=85,
            noise_threshold=0.01,
        )
        similarity = me.calculate_entropy_similarity(spec_query, spec_reference, clean_spectra=False)
        self.assertAlmostEqual(similarity, 1.0, places=4)

    def test_entropy_similarity_ppm(self):
        spec_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)
        spec_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype=np.float32)
        similarity = me.calculate_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_ppm=200)
        self.assertAlmostEqual(similarity, 0.8984397722577456, places=4)

        similarity = me.calculate_entropy_similarity(spec_query, spec_reference, ms2_tolerance_in_ppm=50)
        self.assertAlmostEqual(similarity, 0.8380874603833756, places=4)


if __name__ == "__main__":
    unittest.main()
