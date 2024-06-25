#!/usr/bin/env python3
import json
import numpy as np
from pathlib import Path
from .flash_entropy_search_core import FlashEntropySearchCore


class FlashEntropySearchCoreMediumMemory(FlashEntropySearchCore):
    def __init__(self, path_data, max_ms2_tolerance_in_da=0.024, mz_index_step=0.0001, intensity_weight="entropy") -> None:
        """
        Initialize the EntropySearch class.
        This class is use memmap function to read data from the disk, which is suitable for most of the cases, unless the data is super large.
        
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
