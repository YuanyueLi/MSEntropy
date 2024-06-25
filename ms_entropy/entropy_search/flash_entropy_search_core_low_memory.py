#!/usr/bin/env python3
import json
import numpy as np
from pathlib import Path
from .flash_entropy_search_core import FlashEntropySearchCore


class FlashEntropySearchCoreLowMemory(FlashEntropySearchCore):
    def __init__(self, path_data, max_ms2_tolerance_in_da=0.024, mz_index_step=0.0001, intensity_weight="entropy") -> None:
        """
        Initialize the EntropySearch class.
        This class use file.read function to read the data from the file, which is suitable for very low memory usage.
        
        :param path_data:   The path to save the index data.
        :param max_ms2_tolerance_in_da:  The maximum MS2 tolerance in Da.
        :param mz_index_step:   The step size for the m/z index.
        :param intensity_weight:    The weight for the intensity in the entropy calculation, can be "entropy" or None. Default is "entropy".
            - None: The intensity will not be weighted, then the unweighted similarity will be calculated.
            - "entropy": The intensity will be weighted by the entropy, then the entropy similarity will be calculated.
        """
        super().__init__(max_ms2_tolerance_in_da=max_ms2_tolerance_in_da, mz_index_step=mz_index_step, intensity_weight=intensity_weight)
        self.path_data = Path(str(path_data))
        self.path_data.mkdir(parents=True, exist_ok=True)
        self.index_file = []

    def __del__(self):
        for file in self.index_file:
            file.close()

    def _generate_index_from_peak_data(self, peak_data, max_indexed_mz, append):
        total_peaks_num = peak_data.shape[0]

        # Sort with precursor m/z.
        peak_data.sort(order="ion_mz")

        # Record the m/z, intensity, and spectrum index information for product ions.
        (peak_data["ion_mz"]).tofile(self.path_data / "all_ions_mz.npy")
        (peak_data["intensity"]).tofile(self.path_data / "all_ions_intensity.npy")
        (peak_data["spec_idx"]).tofile(self.path_data / "all_ions_spec_idx.npy")

        # all_ions_mz = self._convert_view_to_array(peak_data.view(np.float32).reshape(total_peaks_num, -1)[:, 0], np.float32, "all_ions_mz")
        # all_ions_intensity = self._convert_view_to_array(peak_data.view(np.float32).reshape(total_peaks_num, -1)[:, 2], np.float32, "all_ions_intensity")
        # all_ions_spec_idx = self._convert_view_to_array(peak_data.view(np.uint32).reshape(total_peaks_num, -1)[:, 3], np.uint32, "all_ions_spec_idx")

        # Assign the index of the product ions.
        peak_data["peak_idx"] = np.arange(0, self.total_peaks_num, dtype=np.uint64)

        # Build index for fast access to the ion's m/z.
        all_ions_mz = np.memmap(self.path_data / "all_ions_mz.npy", dtype=np.float32, mode="r", shape=(total_peaks_num,))
        max_mz = min(np.max(all_ions_mz), max_indexed_mz)
        search_array = np.arange(0.0, max_mz, self.mz_index_step)
        all_ions_mz_idx_start = np.searchsorted(all_ions_mz, search_array, side="left").astype(np.int64)
        all_ions_mz_idx_start.tofile(self.path_data / "all_ions_mz_idx_start.npy")

        ############## Step 3: Build the index by sort with neutral loss mass. ##############
        # Sort with the neutral loss mass.
        peak_data.sort(order="nl_mass")

        # Record the m/z, intensity, spectrum index, and product ions index information for neutral loss ions.
        (peak_data["nl_mass"]).tofile(self.path_data / "all_nl_mass.npy")
        (peak_data["intensity"]).tofile(self.path_data / "all_nl_intensity.npy")
        (peak_data["spec_idx"]).tofile(self.path_data / "all_nl_spec_idx.npy")
        (peak_data["peak_idx"]).tofile(self.path_data / "all_ions_idx_for_nl.npy")

        # all_nl_mass = self._convert_view_to_array(peak_data.view(np.float32).reshape(total_peaks_num, -1)[:, 1], np.float32, "all_nl_mass")
        # all_nl_intensity = self._convert_view_to_array(peak_data.view(np.float32).reshape(total_peaks_num, -1)[:, 2], np.float32, "all_nl_intensity")
        # all_nl_spec_idx = self._convert_view_to_array(peak_data.view(np.uint32).reshape(total_peaks_num, -1)[:, 3], np.uint32, "all_nl_spec_idx")
        # all_ions_idx_for_nl = self._convert_view_to_array(peak_data.view(np.uint64).reshape(total_peaks_num, -1)[:, 2], np.uint64, "all_ions_idx_for_nl")

        # Build the index for fast access to the neutral loss mass.
        all_nl_mass = np.memmap(self.path_data / "all_nl_mass.npy", dtype=np.float32, mode="r", shape=(total_peaks_num,))
        max_mz = min(np.max(all_nl_mass), max_indexed_mz)
        search_array = np.arange(0.0, max_mz, self.mz_index_step)
        all_nl_mass_idx_start = np.searchsorted(all_nl_mass, search_array, side="left").astype(np.int64)
        all_nl_mass_idx_start.tofile(self.path_data / "all_nl_mass_idx_start.npy")

        ############## Step 4: Save the index. ##############
        self.write()
        self.read()
        return self.index

    def read(self, path_data=None):
        """
        Read the index from the file.
        """
        if path_data is not None:
            self.path_data = Path(path_data)

        try:
            self.index = []
            for name in self.index_names:
                self.index.append(np.memmap(self.path_data / f"{name}.npy", dtype=self.index_dtypes[name], mode="r"))
                file_cur = open(self.path_data / f"{name}.npy", "rb")
                file_cur.data_type = np.dtype(self.index_dtypes[name])
                self.index_file.append(file_cur)

            with open(self.path_data / "information.json", "r") as f:
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
        Write the index to the file.
        """
        if path_data is not None:
            assert Path(path_data) == self.path_data, "The path_data is not the same as the path_data in the class."

        information = {
            "mz_index_step": float(self.mz_index_step),
            "total_spectra_num": int(self.total_spectra_num),
            "total_peaks_num": int(self.total_peaks_num),
            "max_ms2_tolerance_in_da": float(self.max_ms2_tolerance_in_da),
        }
        json.dump(information, open(self.path_data / "information.json", "w"))

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
        (
            file_all_ions_mz_idx_start,
            file_all_ions_mz,
            file_all_ions_intensity,
            file_all_ions_spec_idx,
            file_all_nl_mass_idx_start,
            file_all_nl_mass,
            file_all_nl_intensity,
            file_all_nl_spec_idx,
            file_all_ions_idx_for_nl,
        ) = self.index_file
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
                # modified_idx_product = all_ions_spec_idx[product_mz_idx_min:product_mz_idx_max]
                modified_idx_product = _read_data_from_file(file_all_ions_spec_idx, product_mz_idx_min, product_mz_idx_max)

                # modified_value_product = all_ions_intensity[product_mz_idx_min:product_mz_idx_max]
                modified_value_product = _read_data_from_file(file_all_ions_intensity, product_mz_idx_min, product_mz_idx_max)
                modified_value_product = self._score_peaks_with_cpu(intensity, modified_value_product)

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
                # modified_idx_nl = all_nl_spec_idx[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max]
                modified_idx_nl = _read_data_from_file(file_all_nl_spec_idx, neutral_loss_mz_idx_min, neutral_loss_mz_idx_max)
                # modified_value_nl = all_nl_intensity[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max]
                modified_value_nl = _read_data_from_file(file_all_nl_intensity, neutral_loss_mz_idx_min, neutral_loss_mz_idx_max)
                modified_value_nl = self._score_peaks_with_cpu(intensity, modified_value_nl)

                # Check if the neutral loss ion is already matched to other query peak as a product ion
                nl_matched_product_ion_idx = _read_data_from_file(file_all_ions_idx_for_nl, neutral_loss_mz_idx_min, neutral_loss_mz_idx_max)
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
                # modified_idx_product = all_ions_spec_idx[product_mz_idx_min:product_mz_idx_max]
                modified_idx_product = _read_data_from_file(file_all_ions_spec_idx, product_mz_idx_min, product_mz_idx_max)

                # modified_value_product = all_ions_intensity[product_mz_idx_min:product_mz_idx_max]
                modified_value_product = _read_data_from_file(file_all_ions_intensity, product_mz_idx_min, product_mz_idx_max)
                modified_value_product = self._score_peaks_gpu(entropy_transform, intensity, cp.array(modified_value_product))

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
                # print(product_mz_idx_max - product_mz_idx_min, neutral_loss_mz_idx_max - neutral_loss_mz_idx_min)

                # Calculate the entropy similarity for this matched peak
                # modified_idx_nl = all_nl_spec_idx[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max]
                modified_idx_nl = _read_data_from_file(file_all_nl_spec_idx, neutral_loss_mz_idx_min, neutral_loss_mz_idx_max)
                # modified_value_nl = all_nl_intensity[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max]
                modified_value_nl = _read_data_from_file(file_all_nl_intensity, neutral_loss_mz_idx_min, neutral_loss_mz_idx_max)
                modified_value_nl = self._score_peaks_gpu(entropy_transform, intensity, cp.array(modified_value_nl))

                # Check if the neutral loss ion is already matched to other query peak as a product ion
                # nl_matched_product_ion_idx = cp.array(all_ions_idx_for_nl[neutral_loss_mz_idx_min:neutral_loss_mz_idx_max])
                nl_matched_product_ion_idx = cp.array(_read_data_from_file(file_all_ions_idx_for_nl, neutral_loss_mz_idx_min, neutral_loss_mz_idx_max))
                s1 = cp.searchsorted(product_peak_match_idx_min_gpu, nl_matched_product_ion_idx, side="right")
                s2 = cp.searchsorted(product_peak_match_idx_max_gpu, nl_matched_product_ion_idx, side="left")
                modified_value_nl[s1 > s2] = 0

                # Check if this query peak is already matched to a product ion in the same library spectrum
                duplicate_idx_in_nl = self._remove_duplicate_with_gpu(modified_idx_product, modified_idx_nl, self.total_spectra_num)
                modified_value_nl[duplicate_idx_in_nl] = 0

                entropy_similarity_modification_list.append((modified_idx_nl, modified_value_nl.get()))
                del modified_value_nl, nl_matched_product_ion_idx, duplicate_idx_in_nl

            # Release the GPU memory
            cp_dmp = cp.get_default_memory_pool()
            _ = (cp_dmp.used_bytes(), cp_dmp.total_bytes(), cp_dmp.free_bytes())
            cp_dmp.free_all_blocks()
            # print(cp_dmp.used_bytes(), cp_dmp.total_bytes(), cp_dmp.free_bytes(), cp_dmp.get_limit())
            cp.get_default_pinned_memory_pool().free_all_blocks()
            # print(cp_dmp.used_bytes(), cp_dmp.total_bytes(), cp_dmp.free_bytes(), cp_dmp.get_limit())
            entropy_similarity = cp.zeros(self.total_spectra_num, dtype=np.float16)
            for modified_idx, modified_value in entropy_similarity_modification_list:
                entropy_similarity.scatter_add(cp.array(modified_idx), cp.array(modified_value))
            return entropy_similarity.get()
        else:
            raise ValueError("target should be cpu or gpu")


def _read_data_from_file(file_data, item_start, item_end):
    array_type = file_data.data_type
    type_size = array_type.itemsize
    file_data.seek(int(item_start * type_size))
    data = file_data.read(int(type_size * (item_end - item_start)))
    array = np.frombuffer(data, dtype=array_type)
    return array
