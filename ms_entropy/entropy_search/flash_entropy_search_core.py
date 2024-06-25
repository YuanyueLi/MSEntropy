#!/usr/bin/env python3
import json
import numpy as np
from pathlib import Path
from functools import reduce
import multiprocessing
from ..spectra import apply_weight_to_intensity
from .fast_flash_entropy_search import entropy_similarity_search_identity


class FlashEntropySearchCore:
    def __init__(
        self,
        path_data=None,
        max_ms2_tolerance_in_da=0.024,
        mz_index_step=0.0001,
        intensity_weight="entropy", # "entropy" or None
    ) -> None:
        """
        Initialize the EntropySearch class.
    
        :param path_array: The path array of the index files.
        :param max_ms2_tolerance_in_da: The maximum MS2 tolerance used when searching the MS/MS spectra, in Dalton. Default is 0.024.
        :param mz_index_step:   The step size of the m/z index, in Dalton. Default is 0.0001.
                                The smaller the step size, the faster the search, but the larger the index size and longer the index building time.
        :param intensity_weight: The weight of the intensity, can be "entropy" or None. If set to "entropy", the intensity will be weighted by the entropy.
                                If set to None, the intensity will not be weighted, which is equivalent to the unweighted entropy similarity.
        """
        self.mz_index_step = mz_index_step
        self._init_for_multiprocessing = False
        self.max_ms2_tolerance_in_da = max_ms2_tolerance_in_da
        self.intensity_weight = intensity_weight

        self.total_spectra_num = 0
        self.total_peaks_num = 0
        self.index = []

        if path_data:
            self.path_data = Path(path_data)
        else:
            self.path_data = None

        self.index_names = [
            "all_ions_mz_idx_start",
            "all_ions_mz",
            "all_ions_intensity",
            "all_ions_spec_idx",
            "all_nl_mass_idx_start",
            "all_nl_mass",
            "all_nl_intensity",
            "all_nl_spec_idx",
            "all_ions_idx_for_nl",
        ]
        self.index_dtypes = {
            "all_ions_mz_idx_start": np.int64,
            "all_ions_mz": np.float32,
            "all_ions_intensity": np.float32,
            "all_ions_spec_idx": np.uint32,
            "all_nl_mass_idx_start": np.int64,
            "all_nl_mass": np.float32,
            "all_nl_intensity": np.float32,
            "all_nl_spec_idx": np.uint32,
            "all_ions_idx_for_nl": np.uint64,
        }

    def search(
        self,
        method="open",
        target="cpu",
        precursor_mz=None,
        peaks=None,
        ms2_tolerance_in_da=0.02,
        search_type=0,
        search_spectra_idx_min=0,
        search_spectra_idx_max=0,
        output_matched_peak_number=False,
    ):
        """
        Perform identity-, open- or neutral loss search on the MS/MS spectra library.

        :param method:  The search method, can be "open" or "neutral_loss".
                        Set it to "open" for identity search and open search, set it to "neutral_loss" for neutral loss search.
        :param target:  The target to search, can be "cpu" or "gpu".
        :param precursor_mz:    The precursor m/z of the query MS/MS spectrum, required for neutral loss search.
        :param peaks:   The peaks of the query MS/MS spectrum. The peaks need to be precleaned by "clean_spectrum" function.
        :param ms2_tolerance_in_da: The MS2 tolerance used when searching the MS/MS spectra, in Dalton. Default is 0.02.
        :param search_type: The search type, can be 0, 1 or 2.
                            Set it to 0 for searching the whole MS/MS spectra library.
                            Set it to 1 for searching a range of the MS/MS spectra library,
        :param search_spectra_idx_min:  The minimum index of the MS/MS spectra to search, required when search_type is 1.
        :param search_spectra_idx_max:  The maximum index of the MS/MS spectra to search, required when search_type is 1.
        :param output_matched_peak_number: Whether to output the number of matched peaks. Only supported when target is "cpu".
                                            If set to True, the function will return a tuple of (entropy_similarity, matched_peak_number).
        """
        if not self.index:
            return np.zeros(0, dtype=np.float32)
        if len(peaks) == 0:
            return np.zeros(self.total_spectra_num, dtype=np.float32)

        # Check peaks
        assert ms2_tolerance_in_da <= self.max_ms2_tolerance_in_da, "The MS2 tolerance is larger than the maximum MS2 tolerance."
        assert abs(np.sum(peaks[:, 1]) - 1) < 1e-4, "The peaks are not normalized to sum to 1."
        assert (
            peaks.shape[0] <= 1 or np.min(peaks[1:, 0] - peaks[:-1, 0]) > self.max_ms2_tolerance_in_da * 2
        ), "The peaks array should be sorted by m/z, and the m/z difference between two adjacent peaks should be larger than 2 * max_ms2_tolerance_in_da."
        (
            all_ions_mz_idx_start,
            all_ions_mz,
            all_ions_intensity,
            all_ions_spec_idx,
            all_nl_mass_idx_start,
            all_nl_mass,
            all_nl_intensity,
            all_nl_spec_idx,
            all_ions_idx_for_nl,
        ) = self.index
        index_number_in_one_da = int(1 / self.mz_index_step)

        # Prepare the query spectrum
        peaks = self._preprocess_peaks(peaks)

        # Prepare the library
        if method == "open":
            library_mz_idx_start = all_ions_mz_idx_start
            library_mz = all_ions_mz
            library_peaks_intensity = all_ions_intensity
            library_spec_idx = all_ions_spec_idx
        elif method == "neutral_loss":
            library_mz_idx_start = all_nl_mass_idx_start
            library_mz = all_nl_mass
            library_peaks_intensity = all_nl_intensity
            library_spec_idx = all_nl_spec_idx
            peaks[:, 0] = precursor_mz - peaks[:, 0]

        # Start searching
        if target == "cpu":
            entropy_similarity = np.zeros(self.total_spectra_num, dtype=np.float32)
            if output_matched_peak_number:
                matched_peak_number = np.zeros(self.total_spectra_num, dtype=np.uint16)
        else:
            import cupy as cp

            entropy_transform = cp.ElementwiseKernel(
                "T intensity_a, T intensity_b",
                "T similarity",
                """T intensity_ab = intensity_a + intensity_b;
                similarity = intensity_ab * log2f(intensity_ab) - intensity_a * log2f(intensity_a) - intensity_b * log2f(intensity_b);""",
            )
            entropy_similarity = cp.zeros(self.total_spectra_num, dtype=np.float32)

        # Go through all the peaks in the spectrum
        for mz_query, intensity_query in peaks:
            # Determine the mz index range
            product_mz_idx_min = self._find_location_from_array_with_index(
                mz_query - ms2_tolerance_in_da, library_mz, library_mz_idx_start, "left", index_number_in_one_da
            )
            product_mz_idx_max = self._find_location_from_array_with_index(
                mz_query + ms2_tolerance_in_da, library_mz, library_mz_idx_start, "right", index_number_in_one_da
            )

            if target == "cpu" and search_type == 0:
                intensity_library = library_peaks_intensity[product_mz_idx_min:product_mz_idx_max]
                modified_idx = library_spec_idx[product_mz_idx_min:product_mz_idx_max]
                entropy_similarity[modified_idx] += self._score_peaks_with_cpu(intensity_query, intensity_library)
                if output_matched_peak_number:
                    matched_peak_number[modified_idx] += 1
            elif target == "cpu" and search_type == 1:
                entropy_similarity_search_identity(
                    product_mz_idx_min,
                    product_mz_idx_max,
                    intensity_query,
                    entropy_similarity,
                    library_peaks_intensity,
                    library_spec_idx,
                    search_spectra_idx_min,
                    search_spectra_idx_max,
                )
                if output_matched_peak_number:
                    matched_peak_number[library_spec_idx[product_mz_idx_min:product_mz_idx_max]] += 1
            elif target == "gpu":
                intensity_library = cp.array(library_peaks_intensity[product_mz_idx_min:product_mz_idx_max])
                modified_value = entropy_transform(intensity_library, intensity_query)
                modified_idx = cp.array(library_spec_idx[product_mz_idx_min:product_mz_idx_max])
                entropy_similarity.scatter_add(modified_idx, modified_value)

        if target == "cpu":
            if output_matched_peak_number:
                if search_type == 1:
                    matched_peak_number[:search_spectra_idx_min] = 0
                    matched_peak_number[search_spectra_idx_max:] = 0
                return entropy_similarity, matched_peak_number
            else:
                return entropy_similarity
        elif target == "gpu":
            entropy_similarity = entropy_similarity.get()
            if search_type == 1:
                entropy_similarity[:search_spectra_idx_min] = 0
                entropy_similarity[search_spectra_idx_max:] = 0
            return entropy_similarity

    def search_hybrid(self, target="cpu", precursor_mz=None, peaks=None, ms2_tolerance_in_da=0.02):
        """
        Perform the hybrid search for the MS/MS spectra.

        :param target: The target to perform the search. "cpu" for CPU, "gpu" for GPU.
        :param precursor_mz: The precursor m/z of the MS/MS spectra.
        :param peaks: The peaks of the MS/MS spectra, needs to be cleaned with the "clean_spectrum" function.
        :param ms2_tolerance_in_da: The MS/MS tolerance in Da.
        """
        if not self.index:
            return np.zeros(0, dtype=np.float32)
        if len(peaks) == 0:
            return np.zeros(self.total_spectra_num, dtype=np.float32)

        # Check peaks
        assert ms2_tolerance_in_da <= self.max_ms2_tolerance_in_da, "The MS2 tolerance is larger than the maximum MS2 tolerance."
        assert abs(np.sum(peaks[:, 1]) - 1) < 1e-4, "The peaks are not normalized to sum to 1."
        assert (
            peaks.shape[0] <= 1 or np.min(peaks[1:, 0] - peaks[:-1, 0]) > self.max_ms2_tolerance_in_da * 2
        ), "The peaks array should be sorted by m/z, and the m/z difference between two adjacent peaks should be larger than 2 * max_ms2_tolerance_in_da."
        (
            all_ions_mz_idx_start,
            all_ions_mz,
            all_ions_intensity,
            all_ions_spec_idx,
            all_nl_mass_idx_start,
            all_nl_mass,
            all_nl_intensity,
            all_nl_spec_idx,
            all_ions_idx_for_nl,
        ) = self.index
        index_number_in_one_da = int(1 / self.mz_index_step)

        # Prepare the query spectrum
        peaks = self._preprocess_peaks(peaks)

        # Go through all peak in the spectrum and determine the mz index range
        product_peak_match_idx_min = np.zeros(peaks.shape[0], dtype=np.uint64)
        product_peak_match_idx_max = np.zeros(peaks.shape[0], dtype=np.uint64)
        for peak_idx, (mz_query, _) in enumerate(peaks):
            # Determine the mz index range
            product_mz_idx_min = self._find_location_from_array_with_index(
                mz_query - ms2_tolerance_in_da, all_ions_mz, all_ions_mz_idx_start, "left", index_number_in_one_da
            )
            product_mz_idx_max = self._find_location_from_array_with_index(
                mz_query + ms2_tolerance_in_da, all_ions_mz, all_ions_mz_idx_start, "right", index_number_in_one_da
            )

            product_peak_match_idx_min[peak_idx] = product_mz_idx_min
            product_peak_match_idx_max[peak_idx] = product_mz_idx_max

        if target == "cpu":
            entropy_similarity = np.zeros(self.total_spectra_num, dtype=np.float32)
            # Go through all the peaks in the spectrum and calculate the entropy similarity
            for peak_idx, (mz, intensity) in enumerate(peaks):
                ###############################################################
                # Match the original product ion
                product_mz_idx_min = product_peak_match_idx_min[peak_idx]
                product_mz_idx_max = product_peak_match_idx_max[peak_idx]

                # Calculate the entropy similarity for this matched peak
                modified_idx_product = all_ions_spec_idx[product_mz_idx_min:product_mz_idx_max]
                modified_value_product = self._score_peaks_with_cpu(intensity, all_ions_intensity[product_mz_idx_min:product_mz_idx_max])

                entropy_similarity[modified_idx_product] += modified_value_product

                ###############################################################
                # Match the neutral loss ions
                mz_nl = precursor_mz - mz
                # Determine the mz index range
                neutral_loss_mz_idx_min = self._find_location_from_array_with_index(
                    mz_nl - ms2_tolerance_in_da, all_nl_mass, all_nl_mass_idx_start, "left", index_number_in_one_da
                )
                neutral_loss_mz_idx_max = self._find_location_from_array_with_index(
                    mz_nl + ms2_tolerance_in_da, all_nl_mass, all_nl_mass_idx_start, "right", index_number_in_one_da
                )

                # Calculate the entropy similarity for this matched peak
                modified_idx_nl = all_nl_spec_idx[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max]
                modified_value_nl = self._score_peaks_with_cpu(intensity, all_nl_intensity[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max])

                # Check if the neutral loss ion is already matched to other query peak as a product ion
                nl_matched_product_ion_idx = all_ions_idx_for_nl[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max]
                s1 = np.searchsorted(product_peak_match_idx_min, nl_matched_product_ion_idx, side="right")
                s2 = np.searchsorted(product_peak_match_idx_max - 1, nl_matched_product_ion_idx, side="left")
                modified_value_nl[s1 > s2] = 0

                # Check if this query peak is already matched to a product ion in the same library spectrum
                duplicate_idx_in_nl = self._remove_duplicate_with_cpu(modified_idx_product, modified_idx_nl, self.total_spectra_num)
                modified_value_nl[duplicate_idx_in_nl] = 0

                entropy_similarity[modified_idx_nl] += modified_value_nl
            return entropy_similarity

        elif target == "gpu":
            import cupy as cp

            entropy_transform = cp.ElementwiseKernel(
                "T intensity_a, T intensity_b",
                "T similarity",
                """T intensity_ab = intensity_a + intensity_b;
                similarity = intensity_ab * log2f(intensity_ab) - intensity_a * log2f(intensity_a) - intensity_b * log2f(intensity_b);""",
            )
            product_peak_match_idx_min_gpu = cp.array(product_peak_match_idx_min)
            product_peak_match_idx_max_gpu = cp.array(product_peak_match_idx_max - 1)

            entropy_similarity_modification_list = []
            # Go through all the peaks in the spectrum and calculate the entropy similarity
            for peak_idx, (mz, intensity) in enumerate(peaks):
                ###############################################################
                # Match the original product ion
                product_mz_idx_min = product_peak_match_idx_min[peak_idx]
                product_mz_idx_max = product_peak_match_idx_max[peak_idx]

                # Calculate the entropy similarity for this matched peak
                modified_idx_product = all_ions_spec_idx[product_mz_idx_min:product_mz_idx_max]
                modified_value_product = self._score_peaks_gpu(
                    entropy_transform, intensity, cp.array(all_ions_intensity[product_mz_idx_min:product_mz_idx_max])
                )

                entropy_similarity_modification_list.append((modified_idx_product, modified_value_product.get()))
                del modified_value_product

                ###############################################################
                # Match the neutral loss ions
                mz_nl = precursor_mz - mz
                # Determine the mz index range
                neutral_loss_mz_idx_min = self._find_location_from_array_with_index(
                    mz_nl - ms2_tolerance_in_da, all_nl_mass, all_nl_mass_idx_start, "left", index_number_in_one_da
                )
                neutral_loss_mz_idx_max = self._find_location_from_array_with_index(
                    mz_nl + ms2_tolerance_in_da, all_nl_mass, all_nl_mass_idx_start, "right", index_number_in_one_da
                )

                # Calculate the entropy similarity for this matched peak
                modified_idx_nl = all_nl_spec_idx[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max]
                modified_value_nl = self._score_peaks_gpu(
                    entropy_transform, intensity, cp.array(all_nl_intensity[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max])
                )

                # Check if the neutral loss ion is already matched to other query peak as a product ion
                nl_matched_product_ion_idx = cp.array(all_ions_idx_for_nl[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max])
                s1 = cp.searchsorted(product_peak_match_idx_min_gpu, nl_matched_product_ion_idx, side="right")
                s2 = cp.searchsorted(product_peak_match_idx_max_gpu, nl_matched_product_ion_idx, side="left")
                modified_value_nl[s1 > s2] = 0

                # Check if this query peak is already matched to a product ion in the same library spectrum
                duplicate_idx_in_nl = self._remove_duplicate_with_gpu(modified_idx_product, modified_idx_nl, self.total_spectra_num)
                modified_value_nl[duplicate_idx_in_nl] = 0

                entropy_similarity_modification_list.append((modified_idx_nl, modified_value_nl.get()))

            entropy_similarity = cp.zeros(self.total_spectra_num, dtype=np.float32)
            for modified_idx, modified_value in entropy_similarity_modification_list:
                entropy_similarity.scatter_add(cp.array(modified_idx), cp.array(modified_value))
            return entropy_similarity.get()
        else:
            raise ValueError("target should be cpu or gpu")

    def _remove_duplicate_with_cpu(self, array_1, array_2, max_element):
        if len(array_1) + len(array_2) < 4_000_000:
            # When len(array_1) + len(array_2) < 4_000_000, this method is faster than array method
            aux = np.concatenate((array_1, array_2))
            aux_sort_indices = np.argsort(aux, kind="mergesort")
            aux = aux[aux_sort_indices]
            mask = aux[1:] == aux[:-1]
            array2_indices = aux_sort_indices[1:][mask] - array_1.size
            return array2_indices
        else:
            # When len(array_1) + len(array_2) > 4_000_000, this method is faster than sort method
            note = np.zeros(max_element, dtype=np.int8)
            note[array_1] = 1
            duplicate_idx = np.where(note[array_2] == 1)[0]
            return duplicate_idx

    def _remove_duplicate_with_gpu(self, array_1, array_2, max_element):
        import cupy as cp

        if len(array_1) + len(array_2) < 4_000_000:
            # When len(array_1) + len(array_2) < 4_000_000, this method is faster than array method
            aux = cp.array(np.concatenate((array_1, array_2)))
            aux_sort_indices = cp.argsort(aux)
            aux = aux[aux_sort_indices]
            mask = aux[1:] == aux[:-1]
            ar2_indices = aux_sort_indices[1:][mask] - array_1.size
            return ar2_indices

        else:
            # When len(array_1) + len(array_2) > 4_000_000, this method is faster than sort method
            note = cp.zeros(max_element, dtype=np.int8)
            note[array_1] = 1
            duplicate_idx = cp.where(note[array_2] == 1)[0]
            return duplicate_idx

    def build_index(self, all_spectra_list: list, max_indexed_mz: float = 1500.00005, append: bool = False):
        """
        Build the index for the MS/MS spectra library.

        The spectra provided to this function should be a dictionary in the format of {"precursor_mz": precursor_mz, "peaks": peaks}.
        The precursor_mz is the precursor m/z value of the MS/MS spectrum;
        The peaks is a numpy array which has been processed by the function "clean_spectrum".

        :param all_spectra_list:    A list of dictionaries in the format of {"precursor_mz": precursor_mz, "peaks": peaks},
                                    the spectra in the list need to be sorted by the precursor m/z.
        :param max_indexed_mz: The maximum m/z value that will be indexed. Default is 1500.00005.
        :param append:  Not implemented yet.
        """

        # Get the total number of spectra and peaks
        total_peaks_num = int(np.sum([spectrum["peaks"].shape[0] for spectrum in all_spectra_list]))
        total_spectra_num = len(all_spectra_list)
        if append:
            self.total_spectra_num += total_spectra_num
            self.total_peaks_num += total_peaks_num
        else:
            self.total_spectra_num = total_spectra_num
            self.total_peaks_num = total_peaks_num

        # total_spectra_num can not be bigger than 2^32-1 (uint32), total_peak_num can not be bigger than 2^63-1 (int64)
        assert self.total_spectra_num < 2**32 - 1, "The total spectra number is too big."
        assert self.total_peaks_num < 2**63 - 1, "The total peaks number is too big."

        ############## Step 1: Collect the precursor m/z and peaks information. ##############
        peak_data = self._merge_all_spectra_to_peak_data(all_spectra_list, total_peaks_num)

        ############## Step 2: Build the index by sort with product ions. ##############
        self.index = self._generate_index_from_peak_data(peak_data, max_indexed_mz, append=append)
        return self.index

    def _merge_all_spectra_to_peak_data(self, all_spectra_list, total_peaks_num):
        dtype_peak_data = np.dtype(
            [
                ("ion_mz", np.float32),  # The m/z of the fragment ion.
                ("nl_mass", np.float32),  # The neutral loss mass of the fragment ion.
                ("intensity", np.float32),  # The intensity of the fragment ion.
                ("spec_idx", np.uint32),  # The index of the MS/MS spectra.
                ("peak_idx", np.uint64),
            ],
            align=True,
        )  # The index of the fragment ion.

        # Initialize the peak data array.
        peak_data = np.zeros(total_peaks_num, dtype=dtype_peak_data)
        peak_idx = 0

        # Adding the precursor m/z and peaks information to the peak data array.
        for idx, spectrum in enumerate(all_spectra_list):
            precursor_mz, peaks = spectrum["precursor_mz"], spectrum["peaks"]
            # Check the peaks array.
            assert peaks.ndim == 2, "The peaks array should be a 2D numpy array."
            assert peaks.shape[1] == 2, "The peaks array should be a 2D numpy array with the shape of [n, 2]."
            assert peaks.shape[0] > 0, "The peaks array should not be empty."
            assert abs(np.sum(peaks[:, 1]) - 1) < 1e-4, "The peaks array should be normalized to sum to 1."
            assert (
                peaks.shape[0] <= 1 or np.min(peaks[1:, 0] - peaks[:-1, 0]) > self.max_ms2_tolerance_in_da * 2
            ), "The peaks array should be sorted by m/z, and the m/z difference between two adjacent peaks should be larger than 2 * max_ms2_tolerance_in_da."

            # Preprocess the peaks array.
            peaks = self._preprocess_peaks(peaks)

            # Assign the product ion m/z
            peak_data_item = peak_data[peak_idx : (peak_idx + peaks.shape[0])]
            peak_data_item["ion_mz"] = peaks[:, 0]
            # Assign the neutral loss mass
            peak_data_item["nl_mass"] = precursor_mz - peaks[:, 0]
            # Assign the intensity
            peak_data_item["intensity"] = peaks[:, 1]
            # Assign the spectrum index
            peak_data_item["spec_idx"] = idx
            # Set the peak index
            peak_idx += peaks.shape[0]
        return peak_data

    def _generate_index_from_peak_data(self, peak_data, max_indexed_mz, append):
        # Sort with precursor m/z.
        peak_data.sort(order="ion_mz")

        # Record the m/z, intensity, and spectrum index information for product ions.
        all_ions_mz = np.copy(peak_data["ion_mz"])
        all_ions_intensity = np.copy(peak_data["intensity"])
        all_ions_spec_idx = np.copy(peak_data["spec_idx"])

        # Assign the index of the product ions.
        peak_data["peak_idx"] = np.arange(0, self.total_peaks_num, dtype=np.uint64)

        # Build index for fast access to the ion's m/z.
        max_mz = min(np.max(all_ions_mz), max_indexed_mz)
        search_array = np.arange(0.0, max_mz, self.mz_index_step)
        all_ions_mz_idx_start = np.searchsorted(all_ions_mz, search_array, side="left").astype(np.int64)

        ############## Step 3: Build the index by sort with neutral loss mass. ##############
        # Sort with the neutral loss mass.
        peak_data.sort(order="nl_mass")

        # Record the m/z, intensity, spectrum index, and product ions index information for neutral loss ions.
        all_nl_mass = peak_data["nl_mass"]
        all_nl_intensity = peak_data["intensity"]
        all_nl_spec_idx = peak_data["spec_idx"]
        all_ions_idx_for_nl = peak_data["peak_idx"]

        # Build the index for fast access to the neutral loss mass.
        max_mz = min(np.max(all_nl_mass), max_indexed_mz)
        search_array = np.arange(0.0, max_mz, self.mz_index_step)
        all_nl_mass_idx_start = np.searchsorted(all_nl_mass, search_array, side="left").astype(np.int64)

        ############## Step 4: Save the index. ##############
        index = [
            all_ions_mz_idx_start,
            all_ions_mz,
            all_ions_intensity,
            all_ions_spec_idx,
            all_nl_mass_idx_start,
            all_nl_mass,
            all_nl_intensity,
            all_nl_spec_idx,
            all_ions_idx_for_nl,
        ]
        return index

    def _preprocess_peaks(self, peaks):
        """
        Preprocess the peaks.
        """
        if self.intensity_weight == "entropy":
            peaks_clean = np.asarray(apply_weight_to_intensity(peaks))
            peaks_clean[:, 1] /= 2
        elif self.intensity_weight is None:
            peaks_clean = peaks.copy()
            peaks_clean[:, 1] /= 2

        return peaks_clean

    def _score_peaks_with_cpu(self, intensity_query, intensity_library):
        intensity_mix = intensity_library + intensity_query
        modified_value = intensity_mix * np.log2(intensity_mix) - intensity_library * np.log2(intensity_library) - intensity_query * np.log2(intensity_query)
        return modified_value

    def _score_peaks_gpu(self, entropy_transform, intensity_query, intensity_library):
        return entropy_transform(intensity_library, intensity_query)

    def _find_location_from_array_with_index(self, wanted_mz, mz_array, mz_idx_start_array, side, index_number):
        mz_min_int = (np.floor(wanted_mz * index_number)).astype(int)
        mz_max_int = mz_min_int + 1

        if mz_min_int >= len(mz_idx_start_array):
            mz_idx_search_start = mz_idx_start_array[-1]
        else:
            mz_idx_search_start = mz_idx_start_array[mz_min_int].astype(int)

        if mz_max_int >= len(mz_idx_start_array):
            mz_idx_search_end = len(mz_array)
        else:
            mz_idx_search_end = mz_idx_start_array[mz_max_int].astype(int) + 1

        return mz_idx_search_start + np.searchsorted(mz_array[mz_idx_search_start:mz_idx_search_end], wanted_mz, side=side)

    def save_memory_for_multiprocessing(self):
        """
        Move the numpy array in the index to shared memory in order to save memory.
        This function is not required when you only use one thread to search the MS/MS spectra.
        When use multiple threads, this function is also not required but highly recommended, as it avoids the memory copy and saves a lot of memory and time.
        """
        if self._init_for_multiprocessing:
            return

        for i, array in enumerate(self.index):
            self.index[i] = _convert_numpy_array_to_shared_memory(array)
        self._init_for_multiprocessing = True

    def read(self, path_data=None):
        """
        Read the index from the specified path.
        """
        try:
            if path_data is None:
                path_data = self.path_data

            path_data = Path(path_data)
            self.index = []
            for name in self.index_names:
                self.index.append(np.fromfile(path_data / f"{name}.npy", dtype=self.index_dtypes[name]))

            with open(path_data / "information.json", "r") as f:
                information = json.load(f)
            self.mz_index_step = information["mz_index_step"]
            self.total_spectra_num = information["total_spectra_num"]
            self.total_peaks_num = information["total_peaks_num"]
            self.max_ms2_tolerance_in_da = information["max_ms2_tolerance_in_da"]
            return True
        except:
            return False

    def write(self, path_data=None):
        """
        Write the index to the specified path.
        """
        if path_data is None:
            path_data = self.path_data

        path_data = Path(path_data)
        path_data.mkdir(parents=True, exist_ok=True)
        for i, name in enumerate(self.index_names):
            self.index[i].tofile(str(path_data / f"{name}.npy"))
        information = {
            "mz_index_step": float(self.mz_index_step),
            "total_spectra_num": int(self.total_spectra_num),
            "total_peaks_num": int(self.total_peaks_num),
            "max_ms2_tolerance_in_da": float(self.max_ms2_tolerance_in_da),
        }
        with open(path_data / "information.json", "w") as f:
            json.dump(information, f)


def _convert_numpy_array_to_shared_memory(np_array, array_c_type=None):
    """
    The char table of shared memory can be find at:
    https://docs.python.org/3/library/struct.html#format-characters
    https://docs.python.org/3/library/array.html#module-array (This one is wrong!)
    The documentation of numpy.frombuffer can be find at:
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.frombuffer.html
    Note: the char table is different from the char table in numpy
    """
    dim = np_array.shape
    num = reduce(lambda x, y: x * y, dim)
    if array_c_type is None:
        array_c_type = np_array.dtype.char
    base = multiprocessing.Array(array_c_type, num, lock=False)
    np_array_new = np.frombuffer(base, dtype=np_array.dtype).reshape(dim)
    np_array_new[:] = np_array
    return np_array_new
