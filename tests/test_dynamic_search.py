import numpy as np
import unittest
import tempfile
from ms_entropy import FlashEntropySearchCoreForDynamicIndexing, clean_spectrum


class TestFlashEntropySearchWithCpu(unittest.TestCase):
    def setUp(self):
        path_test = tempfile.mkdtemp()
        spectral_library = [
            {"id": "Demo spectrum 1", "precursor_mz": 150.0, "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0], [103.0, 1.0]], dtype=np.float32)},
            {
                "id": "Demo spectrum 2",
                "precursor_mz": 220.0,
                "peaks": np.array([[200.0, 1.0], [101.0, 1.0], [202.0, 1.0], [204.0, 1.0], [205.0, 1.0]], dtype=np.float32),
            },
            {
                "id": "Demo spectrum 3",
                "precursor_mz": 250.0,
                "peaks": np.array([[100.0, 1.0], [201.0, 1.0], [202.0, 1.0], [104.0, 1.0], [105.0, 1.0]], dtype=np.float32),
            },
            {"id": "Demo spectrum 4", "precursor_mz": 350.0, "peaks": [[100.0, 1.0], [101.0, 1.0], [302.0, 1.0], [104.0, 1.0], [105.0, 1.0]]},
        ]
        query_spectrum = {"precursor_mz": 150.0, "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0], [103.0, 1.0]], dtype=np.float32)}

        for spec in spectral_library:
            # Clean the peaks
            spec["peaks"] = clean_spectrum(
                peaks=spec["peaks"],
                max_mz=spec["precursor_mz"],
                noise_threshold=0.01,
                min_ms2_difference_in_da=0.02,
            )

        self.flash_entropy = FlashEntropySearchCoreForDynamicIndexing(path_data=path_test)
        self.flash_entropy.build_index(spectral_library)
        query_spectrum["peaks"] = clean_spectrum(max_mz=query_spectrum["precursor_mz"], peaks=query_spectrum["peaks"])
        self.query_spectrum = query_spectrum
        self.spectral_library = spectral_library

    def test_read_and_write(self):
        self.flash_entropy.write()
        self.flash_entropy.read()

    def test_hybrid_search(self):
        similarity = self.flash_entropy.search_hybrid(
            precursor_mz=self.query_spectrum["precursor_mz"], peaks=self.query_spectrum["peaks"], ms2_tolerance_in_da=0.02
        )
        np.testing.assert_almost_equal(similarity, [1.0, 0.22299, 0.66897, 0.66897], decimal=5)

    def test_neutral_loss_search(self):
        similarity = self.flash_entropy.search(
            precursor_mz=self.query_spectrum["precursor_mz"], peaks=self.query_spectrum["peaks"], method="neutral_loss", ms2_tolerance_in_da=0.02
        )
        np.testing.assert_almost_equal(similarity, [1.0, 0.0, 0.44598, 0.22299], decimal=5)
        similarity, matched_peaks = self.flash_entropy.search(
            precursor_mz=self.query_spectrum["precursor_mz"],
            peaks=self.query_spectrum["peaks"],
            method="neutral_loss",
            ms2_tolerance_in_da=0.02,
            output_matched_peak_number=True,
        )
        np.testing.assert_almost_equal(similarity, [1.0, 0.0, 0.44598, 0.22299], decimal=5)
        np.testing.assert_almost_equal(matched_peaks, [4, 0, 2, 1], decimal=5)

    def test_open_search(self):
        similarity = self.flash_entropy.search(peaks=self.query_spectrum["peaks"], method="open", ms2_tolerance_in_da=0.02)
        np.testing.assert_almost_equal(similarity, [1.0, 0.22299, 0.22299, 0.44598], decimal=5)
        similarity, matched_peaks = self.flash_entropy.search(
            peaks=self.query_spectrum["peaks"], method="open", ms2_tolerance_in_da=0.02, output_matched_peak_number=True
        )
        np.testing.assert_almost_equal(similarity, [1.0, 0.22299, 0.22299, 0.44598], decimal=5)
        np.testing.assert_almost_equal(matched_peaks, [4, 1, 1, 2], decimal=5)

    def test_identity_search(self):
        similarity = self.flash_entropy.search(method="open", peaks=self.query_spectrum["peaks"], ms2_tolerance_in_da=0.02)
        library_precursor_mz = np.array([x["precursor_mz"] for x in self.spectral_library], dtype=np.float32)
        similarity[np.abs(library_precursor_mz - self.query_spectrum["precursor_mz"]) > 0.01] = 0
        np.testing.assert_almost_equal(similarity, [1.0, 0.0, 0.0, 0.0], decimal=5)
        similarity, matched_peaks = self.flash_entropy.search(
            method="open",
            peaks=self.query_spectrum["peaks"],
            ms2_tolerance_in_da=0.02,
            output_matched_peak_number=True,
        )
        similarity[np.abs(library_precursor_mz - self.query_spectrum["precursor_mz"]) > 0.01] = 0
        matched_peaks[np.abs(library_precursor_mz - self.query_spectrum["precursor_mz"]) > 0.01] = 0
        np.testing.assert_almost_equal(similarity, [1.0, 0.0, 0.0, 0.0], decimal=5)
        np.testing.assert_almost_equal(matched_peaks, [4, 0, 0, 0], decimal=5)


if __name__ == "__main__":
    unittest.main()
