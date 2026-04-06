import numpy as np
from ..spectra import apply_weight_to_intensity
from pathlib import Path
import json


class DynamicEntropySearchCore(object):
    def __init__(
        self,
        path_data,
        max_ms2_tolerance_in_da=0.024,
        extend_fold=3,
        mass_per_block: float = 0.05,
        max_indexed_mz: float = 1500.00005,
        intensity_weight="entropy",  # "entropy" or None
    ) -> None:
        
        
        """
        Initialize the :class:`DynamicEntropySearchCore` object.

        Parameters
        ----------
        path_data : str or Path
            Path to the directory where index data for this group will be stored.

        max_ms2_tolerance_in_da : float, optional
            Maximum MS2 fragment-ion tolerance (in Daltons) allowed during entropy search. 
            Default is ``0.024``.

        extend_fold : int, optional
            Expansion factor applied when allocating space for block data.
            Must be greater than 1. 
            Default is ``3``.

        mass_per_block : float, optional
            The m/z interval used to define ion index blocks.  
            Default is ``0.05`` Da.

        max_indexed_mz : float, optional
            Maximum fragment-ion m/z to index. Fragment ions with larger m/z values are placed into a single block.  
            Default is ``1500.00005``.

        intensity_weight : {"entropy", None}, optional
            Whether fragment intensities should be entropy-weighted.  
            If ``"entropy"``, intensities are transformed using entropy weighting.  
            If ``None``, raw intensity values are used.  
            Default is ``"entropy"``.


        Returns
        -------
        None
        """
        self.max_ms2_tolerance_in_da = max_ms2_tolerance_in_da
        self.mass_per_block = mass_per_block
        self.max_indexed_mz = max_indexed_mz

        assert extend_fold > 1, "The extend_fold should be larger than 1."
        self.extend_fold = extend_fold
        self.intensity_weight = intensity_weight

        self.total_spectra_num = 0
        self.total_peaks_num = 0

        self.path_data = Path(path_data)
        self.path_data.mkdir(parents=True, exist_ok=True)

        self.index = []

        self.dtype_block_data = np.dtype(
            [
                ("spec_idx", np.uint32),  # The index of the MS/MS spectra.
                ("mass", np.float32),  # The m/z of the fragment ion.
                ("intensity", np.float32),  # The intensity of the fragment ion.
            ],
            align=True,
        )  # The index of the fragment ion.

        self.dtype_block_data_nl = np.dtype(
            [
                ("fragment_mz", np.float32),  # The m/z of corresponding ion
                ("spec_idx", np.uint32),  # The index of the MS/MS spectra.
                ("mass", np.float32),  # The mass of the neutral loss.
                ("intensity", np.float32),  # The intensity of the neutral loss.
            ],
            align=True,
        )  # The index of the neutral loss.

        self.dtype_block_info = np.dtype(
            [
                ("start_idx", np.uint64),  # The start index of the block.
                ("data_len", np.uint32),  # The length of data in the block.
                ("reserved_len", np.uint32),  # The reserved length of the block.
                ("is_sorted", bool),  # Whether the block is sorted or not
            ]
        )
        self.index_names = [
            "block_ions_info",
            "block_nl_info",
        ]  # The information of block_data and block_data_nl

        self.index_dtypes = {
            "block_ions_info": self.dtype_block_info,
            "block_nl_info": self.dtype_block_info,
        }

    def search(
        self,
        method="open",
        precursor_mz=None,
        peaks=None,
        ms2_tolerance_in_da=0.02,
    ):
        
        """
        Perform open search or neutral-loss search on the MS/MS spectral library.

        Parameters
        ----------
        method : {"open", "neutral_loss"}, optional
            Search mode.

            - ``"open"`` — identity search or open search  
            - ``"neutral_loss"`` — neutral-loss-based matching

            Default is ``"open"``.

        precursor_mz : float, optional
            Precursor m/z of the query MS/MS spectrum.  
            Required when ``method="neutral_loss"``.

        peaks : numpy.ndarray
            Array of fragment peaks from the query spectrum.  
            The peaks must be preprocessed using :func:`clean_spectrum` and normalized such that intensities sum to 1.

        ms2_tolerance_in_da : float, optional
            Fragment mass tolerance (Da) used for peak matching.  
            Must be less than or equal to ``max_ms2_tolerance_in_da``.  
            Default is ``0.02``.

        Returns
        -------
        numpy.ndarray
            A 1-D array of entropy similarity scores for all spectra in the library, with dtype ``float32``. 
            Each entry corresponds to the similarity score for one reference spectrum.


        """

        if not self.index or len(self.index) != 2:
            raise RuntimeError("Index not loaded. Call build_index(...) or read(...) first.")
        if len(peaks) == 0:
            return np.zeros(self.total_spectra_num, dtype=np.float32)

        # Check peaks
        assert ms2_tolerance_in_da <= self.max_ms2_tolerance_in_da, "The MS2 tolerance is larger than the maximum MS2 tolerance."
        assert abs(np.sum(peaks[:, 1]) - 1) < 1e-4, "The peaks are not normalized to sum to 1."
        assert (
            peaks.shape[0] <= 1 or np.min(peaks[1:, 0] - peaks[:-1, 0]) > self.max_ms2_tolerance_in_da * 2
        ), "The peaks array should be sorted by m/z, and the m/z difference between two adjacent peaks should be larger than 2 * max_ms2_tolerance_in_da."

        # Prepare the query spectrum
        peaks=np.asarray(peaks, dtype=np.float32)
        peaks = self._preprocess_peaks(peaks)

        [block_product_info, block_nl_info] = self.index

        if method == "neutral_loss" and block_nl_info is None:
            raise RuntimeError("Neutral loss search cannot be implemented because neutral loss index is not built.")

        # Prepare the library
        if method == "open":
            # open search library
            block_ions_data = np.memmap(self.path_data / "ions_data.bin", dtype=self.dtype_block_data, mode="r")
            block_ions_info = block_product_info

        elif method == "neutral_loss":
            if precursor_mz is None:
                raise ValueError("precursor_mz should not be None in 'neutral_loss' method.")
            # neutral loss search library
            block_ions_data = np.memmap(self.path_data / "nl_data.bin", dtype=self.dtype_block_data_nl, mode="r")
            block_ions_info = block_nl_info
            # generate the nl_mass of the query spectrum
            peaks[:, 0] = float(precursor_mz) - peaks[:, 0]
        else:
            raise ValueError(f"Unknown method: {method}. Use 'open' or 'neutral_loss'.")

        # Start searching
        entropy_similarity = np.zeros(self.total_spectra_num, dtype=np.float32)

        # find the block location of query peaks; same process for both open search and neutral loss search
        min_block_query_idx, max_block_query_idx = self._locate_query_peaks(peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)
        # Go through all the peaks in the spectrum
        for i, (mz_query, intensity_query) in enumerate(peaks):
            # Determine the mz index range
            wanted_block_left = min_block_query_idx[i]
            wanted_block_right = max_block_query_idx[i]

            for block in range(wanted_block_left, (wanted_block_right + 1)):
                block_start_idx = block_ions_info["start_idx"][block]
                block_length = block_ions_info["data_len"][block]
                block_is_sorted = block_ions_info["is_sorted"][block]
                selected_block_data = block_ions_data[block_start_idx : (block_start_idx + block_length)]
                selected_block_data_mass = selected_block_data["mass"]
                if block_is_sorted == False:
                    select_filter = np.bitwise_and(
                        (mz_query - ms2_tolerance_in_da <= selected_block_data_mass), (selected_block_data_mass <= mz_query + ms2_tolerance_in_da)
                    )
                    block_spec_idx = selected_block_data["spec_idx"][select_filter]
                    block_intensity = selected_block_data["intensity"][select_filter]

                    if np.all(select_filter == 0):
                        continue

                elif block_is_sorted == True:
                    left_idx = np.searchsorted(selected_block_data_mass, mz_query - ms2_tolerance_in_da, side="left")
                    right_idx = np.searchsorted(selected_block_data_mass, mz_query + ms2_tolerance_in_da, side="right")

                    if left_idx == right_idx:
                        continue
                    block_spec_idx = selected_block_data["spec_idx"][left_idx:right_idx]
                    block_intensity = selected_block_data["intensity"][left_idx:right_idx]

                entropy_similarity[block_spec_idx] += self._score_peaks_with_cpu(intensity_query, block_intensity)

        return entropy_similarity

    def search_hybrid(self, precursor_mz, peaks, ms2_tolerance_in_da=0.02):

        """
        Perform hybrid search against the MS/MS spectral index.

        Hybrid search incorporates both:
        - **Open search** (direct fragment ion matching)
        - **Neutral-loss search** (matching peaks transformed as ``precursor_mz - peak_mz``)

        Parameters
        ----------
        precursor_mz : float
            The precursor m/z of the query MS/MS spectrum.

        peaks : numpy.ndarray
            Fragment ions of the query spectrum.  
            Must be preprocessed using :func:`clean_spectrum`, sorted by m/z, and normalized such that the sum of intensities equals 1.
        
        ms2_tolerance_in_da : float, optional
            Mass tolerance (Da) for fragment matching.  
            Must be ≤ ``max_ms2_tolerance_in_da``.  
            Default is ``0.02``.

        Returns
        -------
        numpy.ndarray
            A 1-D array of entropy similarity scores between the query spectrum and each library spectrum.  
            Dtype is ``float32`` with length equal to the number of indexed spectra.


        """

        if not self.index:
            raise RuntimeError("Index not loaded. Call build_index(...) or read(...) first.")

        if len(peaks) == 0:
            return np.zeros(self.total_spectra_num, dtype=np.float32)

        # Check peaks
        assert ms2_tolerance_in_da <= self.max_ms2_tolerance_in_da, "The MS2 tolerance is larger than the maximum MS2 tolerance."
        assert abs(np.sum(peaks[:, 1]) - 1) < 1e-4, "The peaks are not normalized to sum to 1."
        assert (
            peaks.shape[0] <= 1 or np.min(peaks[1:, 0] - peaks[:-1, 0]) > self.max_ms2_tolerance_in_da * 2
        ), "The peaks array should be sorted by m/z, and the m/z difference between two adjacent peaks should be larger than 2 * max_ms2_tolerance_in_da."

        # Prepare the query spectrum
        peaks=np.asarray(peaks, dtype=np.float32)
        peaks = self._preprocess_peaks(peaks)

        [block_ions_info, block_nl_info] = self.index

        ###############################################################
        # Match the original product ion
        # find the block location of query peaks for open search
        min_block_query_idx, max_block_query_idx = self._locate_query_peaks(peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)

        block_ions_data = np.memmap(self.path_data / "ions_data.bin", dtype=self.dtype_block_data, mode="r")

        ###############################################################
        # Match the neutral loss ions
        # find the block location of query peaks for neutral loss search
        nl_peaks = np.copy(peaks)
        nl_peaks[:, 0] = float(precursor_mz) - peaks[:, 0]
        nl_min_block_query_idx, nl_max_block_query_idx = self._locate_query_peaks(peaks=nl_peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)

        block_nl_data = np.memmap(self.path_data / "nl_data.bin", dtype=self.dtype_block_data_nl, mode="r")

        # collect the ions information based on the block_query_idx

        entropy_similarity = np.zeros(self.total_spectra_num, dtype=np.float32)

        all_ions_mz_left = np.array(peaks[:, 0] - ms2_tolerance_in_da, dtype=np.float32)
        all_ions_mz_right = np.array(peaks[:, 0] + ms2_tolerance_in_da, dtype=np.float32)

        for i, (mz_query, intensity_query) in enumerate(peaks):
            wanted_block_left = min_block_query_idx[i]
            wanted_block_right = max_block_query_idx[i]

            matched_spec_idx = []
            for block in range(wanted_block_left, (wanted_block_right + 1)):

                block_start_idx = block_ions_info["start_idx"][block]
                block_length = block_ions_info["data_len"][block]
                block_is_sorted = block_ions_info["is_sorted"][block]
                selected_block_data = block_ions_data[block_start_idx : (block_start_idx + block_length)]
                selected_block_data_mass = selected_block_data["mass"]

                if block_is_sorted == False:
                    select_filter = np.bitwise_and(
                        (mz_query - ms2_tolerance_in_da <= selected_block_data_mass), (selected_block_data_mass <= mz_query + ms2_tolerance_in_da)
                    )

                    if np.all(select_filter == 0):
                        continue

                    block_spec_idx = selected_block_data["spec_idx"][select_filter]
                    block_intensity = selected_block_data["intensity"][select_filter]

                elif block_is_sorted == True:

                    left_idx = np.searchsorted(selected_block_data_mass, mz_query - ms2_tolerance_in_da, side="left")
                    right_idx = np.searchsorted(selected_block_data_mass, mz_query + ms2_tolerance_in_da, side="right")

                    if left_idx == right_idx:
                        continue
                    block_spec_idx = selected_block_data["spec_idx"][left_idx:right_idx]
                    block_intensity = selected_block_data["intensity"][left_idx:right_idx]
                entropy_similarity[block_spec_idx] += self._score_peaks_with_cpu(intensity_query, block_intensity)

                matched_spec_idx.append(block_spec_idx)

            if matched_spec_idx:
                matched_spectra_idx = np.concatenate(matched_spec_idx, dtype=np.uint32)
            else:
                matched_spectra_idx = np.zeros(shape=(0,), dtype=np.uint32)

            nl_mass = float(precursor_mz) - mz_query
            nl_wanted_block_left = nl_min_block_query_idx[i]
            nl_wanted_block_right = nl_max_block_query_idx[i]

            nl_spec_idx = []
            nl_intensity = []
            nl_fragment_mz = []
            for nl_block in range(nl_wanted_block_left, (nl_wanted_block_right + 1)):
                nl_block_start_idx = block_nl_info["start_idx"][nl_block]
                nl_block_length = block_nl_info["data_len"][nl_block]
                nl_block_is_sorted = block_nl_info["is_sorted"][nl_block]
                selected_nl_data = block_nl_data[nl_block_start_idx : (nl_block_start_idx + nl_block_length)]
                selected_nl_data_mass = selected_nl_data["mass"]

                if nl_block_is_sorted == False:
                    select_nl_filter = np.bitwise_and(
                        (nl_mass - ms2_tolerance_in_da <= selected_nl_data_mass), (selected_nl_data_mass <= nl_mass + ms2_tolerance_in_da)
                    )

                    if np.all(select_nl_filter == 0):
                        continue

                    nl_block_spec_idx = selected_nl_data["spec_idx"][select_nl_filter]
                    nl_block_intensity = selected_nl_data["intensity"][select_nl_filter]
                    nl_block_fragment_mz = selected_nl_data["fragment_mz"][select_nl_filter]

                elif nl_block_is_sorted == True:

                    left_nl_idx = np.searchsorted(selected_nl_data_mass, nl_mass - ms2_tolerance_in_da, side="left")
                    right_nl_idx = np.searchsorted(selected_nl_data_mass, nl_mass + ms2_tolerance_in_da, side="right")

                    if left_nl_idx == right_nl_idx:
                        continue

                    nl_block_spec_idx = selected_nl_data["spec_idx"][left_nl_idx:right_nl_idx]
                    nl_block_intensity = selected_nl_data["intensity"][left_nl_idx:right_nl_idx]
                    nl_block_fragment_mz = selected_nl_data["fragment_mz"][left_nl_idx:right_nl_idx]

                nl_spec_idx.append(nl_block_spec_idx)
                nl_intensity.append(nl_block_intensity)
                nl_fragment_mz.append(nl_block_fragment_mz)

            if nl_spec_idx:
                nl_all_spectra_idx = np.concatenate(nl_spec_idx)
                nl_all_fragment_mz = np.concatenate(nl_fragment_mz)
                nl_all_intensity = np.concatenate(nl_intensity)

                # Check if the neutral loss ion is already matched to other query peak as a product ion
                modified_value_nl = self._score_peaks_with_cpu(intensity_query, nl_all_intensity)
                f1 = np.searchsorted(all_ions_mz_left, nl_all_fragment_mz, side="right")
                f2 = np.searchsorted(all_ions_mz_right, nl_all_fragment_mz, side="left")
                modified_value_nl[f1 > f2] = 0

                # Check if this query peak is already matched to a product ion in the same library spectrum
                duplicate_idx_in_nl = self._remove_duplicate_with_cpu(matched_spectra_idx, nl_all_spectra_idx, self.total_spectra_num)

                modified_value_nl[duplicate_idx_in_nl] = 0
                entropy_similarity[nl_all_spectra_idx] += modified_value_nl

        return entropy_similarity

    def _locate_query_peaks(self, peaks: np.ndarray, ms2_tolerance_in_da: np.float32):
        # Find the block location of query peaks

        search_array = np.arange(0.0, self.max_indexed_mz, self.mass_per_block)
        blocks_num = len(search_array)
        last_idx = blocks_num - 1

        min_block_query_idx = np.zeros(peaks.shape[0], dtype=np.int64)
        max_block_query_idx = np.zeros(peaks.shape[0], dtype=np.int64)

        for peak_idx, (mz_query, _) in enumerate(peaks):
            left = int((mz_query - ms2_tolerance_in_da) / self.mass_per_block)
            right = int((mz_query + ms2_tolerance_in_da) / self.mass_per_block)

            left = max(0, min(left, last_idx))
            right = max(0, min(right, last_idx))

            min_block_query_idx[peak_idx] = left
            max_block_query_idx[peak_idx] = right

        return min_block_query_idx, max_block_query_idx

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

    def build_index(self, all_spectra_list: list, index_for_neutral_loss: bool = True):

        """
        Build the fragment-ion index (and optional neutral-loss index) for the MS/MS spectral library.

        The input spectra must be dictionaries with the structure::

            {
                "precursor_mz": float,
                "peaks": numpy.ndarray
            }

        where:

        - ``precursor_mz`` is the precursor m/z value of the MS/MS spectrum  
        - ``peaks`` is a 2-column numpy array that has been preprocessed using :func:`clean_spectrum`, containing sorted and normalized (m/z, intensity) pairs

        Parameters
        ----------
        all_spectra_list : list of dict
            A list of spectrum dictionaries in the format described above.  
            
        index_for_neutral_loss : bool, optional
            Whether to also build the neutral-loss index.  
            If ``True``, the method generates both the product-ion and neutral-loss indices.  
            If ``False``, only the product-ion index is built.  
            Default is ``True``.

        Returns
        -------
        None
            The method updates internal index structures in-place.

        """

        # Get the total number of spectra and peaks
        total_peaks_num = int(np.sum([np.array(spectrum["peaks"]).shape[0] for spectrum in all_spectra_list]))
        total_spectra_num = len(all_spectra_list)

        self.total_spectra_num = total_spectra_num
        self.total_peaks_num = total_peaks_num

        # total_spectra_num can not be bigger than 2^32-1 (uint32), total_peak_num can not be bigger than 2^63-1 (int64)
        assert self.total_spectra_num < 2**32 - 1, "The total spectra number is too big."
        assert self.total_peaks_num < 2**63 - 1, "The total peaks number is too big."

        ############## Step 1: Collect the precursor m/z and peaks information. ##############
        peak_data = self._merge_all_spectra_to_peak_data(all_spectra_list, total_peaks_num)

        ############## Step 2: Build the index by sort with product ions. ##############
        self._generate_index_from_peak_data(peak_data, index_for_neutral_loss=index_for_neutral_loss)

        return

    def _merge_all_spectra_to_peak_data(
        self,
        all_spectra_list,
        total_peaks_num,
        append: bool = False,  # when it is a new spectra_list needed to be added into the old index, append is True
        existed_max_spec_idx: int = None,  # when it is a new spectra_list needed to be added into the old index, existed_max_spec_idx is the max num of existing spectra
    ):

        dtype_peak_data = np.dtype(
            [
                ("ion_mz", np.float32),  # The m/z of the fragment ion.
                ("nl_mass", np.float32),  # The neutral loss mass of the fragment ion.
                ("intensity", np.float32),  # The intensity of the fragment ion.
                ("spec_idx", np.uint32),  # The index of the MS/MS spectra.
            ],
            align=True,
        )  # The index of the fragment ion.

        # Initialize the peak data array.
        peak_data = np.zeros(total_peaks_num, dtype=dtype_peak_data)
        peak_idx = 0
        # Adding the precursor m/z and peaks information to the peak data array.
        for idx, spectrum in enumerate(all_spectra_list):
            precursor_mz, peaks = spectrum["precursor_mz"], np.array(spectrum["peaks"])
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
            peak_data_item["nl_mass"] = float(precursor_mz) - peaks[:, 0]
            # Assign the intensity
            peak_data_item["intensity"] = peaks[:, 1]
            # Assign the spectrum index
            if append:
                peak_data_item["spec_idx"] = existed_max_spec_idx + 1 + idx
                self.total_spectra_num += 1
                self.total_peaks_num += peaks.shape[0]
            else:
                peak_data_item["spec_idx"] = idx
            # Set the peak index
            peak_idx += peaks.shape[0]

        return peak_data

    def _generate_index_from_peak_data(self, peak_data, index_for_neutral_loss: bool = True):

        extend_fold = self.extend_fold
        ### Open Search ###
        # Sort with precursor m/z
        peak_data.sort(order="ion_mz")
        # Assign the index of the product ions.

        search_array = np.arange(0.0, self.max_indexed_mz, self.mass_per_block)
        # Calculate the number of blocks
        blocks_num = len(search_array)

        with open(self.path_data / "ions_data.bin", "wb") as f_block_ions_data:
            # Assign the open search block
            block_ions_info = np.empty(blocks_num, dtype=self.dtype_block_info)
            block_ions_info["is_sorted"] = True * blocks_num

            block_data_cur_idx = 0
            block_start_loc_array = np.searchsorted(peak_data["ion_mz"], search_array, side="left")

            for block_idx in range(blocks_num):
                block_start_loc = block_start_loc_array[block_idx]
                if block_idx < blocks_num - 1:
                    block_end_loc = block_start_loc_array[block_idx + 1]
                else:
                    block_end_loc = len(peak_data)

                # Write block information
                data_len = block_end_loc - block_start_loc
                reserved_len = int(data_len * extend_fold)
                block_ions_info[block_idx]["start_idx"] = block_data_cur_idx
                block_ions_info[block_idx]["data_len"] = data_len
                block_ions_info[block_idx]["reserved_len"] = reserved_len

                if block_end_loc == block_start_loc:
                    # If there are no peaks in this block, skip it.
                    continue

                # Generate block data
                block_data = np.empty(reserved_len, dtype=self.dtype_block_data)
                block_data["spec_idx"][:data_len] = peak_data["spec_idx"][block_start_loc:block_end_loc]
                block_data["mass"][:data_len] = peak_data["ion_mz"][block_start_loc:block_end_loc]
                block_data["intensity"][:data_len] = peak_data["intensity"][block_start_loc:block_end_loc]

                # Write block data

                block_data.tofile(f_block_ions_data)

                # Move to the next block
                block_data_cur_idx += reserved_len

        ### Neutral Loss Search ###
        if index_for_neutral_loss:
            # Sort with the neutral loss mass.
            peak_data.sort(order="nl_mass")
            # Record the nl_mass, intensity, spectrum index, and product ions index information for neutral loss ions.
            with open(self.path_data / "nl_data.bin", "wb") as f_block_nl_data:
                # Assign the neutral_loss search block
                block_nl_info = np.empty(blocks_num, dtype=self.dtype_block_info)

                block_nl_info["is_sorted"] = True * blocks_num

                block_data_cur_idx = 0
                block_start_loc_array = np.searchsorted(peak_data["nl_mass"], search_array, side="left")

                for block_idx in range(blocks_num):
                    block_start_loc = block_start_loc_array[block_idx]
                    if block_idx < blocks_num - 1:
                        block_end_loc = block_start_loc_array[block_idx + 1]
                    else:
                        block_end_loc = len(peak_data)

                    # Write block information
                    data_len = block_end_loc - block_start_loc
                    reserved_len = int(data_len * extend_fold)
                    block_nl_info[block_idx]["start_idx"] = block_data_cur_idx
                    block_nl_info[block_idx]["data_len"] = data_len
                    block_nl_info[block_idx]["reserved_len"] = reserved_len

                    if block_end_loc == block_start_loc:
                        # If there are no peaks in this block, skip it.
                        continue

                    # Generate block data
                    block_data = np.empty(reserved_len, dtype=self.dtype_block_data_nl)
                    block_data["spec_idx"][:data_len] = peak_data["spec_idx"][block_start_loc:block_end_loc]
                    block_data["mass"][:data_len] = peak_data["nl_mass"][block_start_loc:block_end_loc]
                    block_data["intensity"][:data_len] = peak_data["intensity"][block_start_loc:block_end_loc]
                    block_data["fragment_mz"][:data_len] = peak_data["ion_mz"][block_start_loc:block_end_loc]

                    # Write block data

                    block_data.tofile(f_block_nl_data)

                    # Move to the next block
                    block_data_cur_idx += reserved_len

        else:
            block_nl_info = None

        self.index = [
            block_ions_info,
            block_nl_info,
        ]

        return

    def _preprocess_peaks(self, peaks):
        """
        Preprocess the peaks.
        """
        if self.intensity_weight == "entropy":
            peaks_clean = np.asarray(apply_weight_to_intensity(peaks))
        elif self.intensity_weight is None:
            peaks_clean = peaks.copy()
        peaks_clean[:, 1] /= 2

        return peaks_clean

    def _score_peaks_with_cpu(self, intensity_query, intensity_library):
        intensity_mix = intensity_library + intensity_query
        modified_value = intensity_mix * np.log2(intensity_mix) - intensity_library * np.log2(intensity_library) - intensity_query * np.log2(intensity_query)
        return modified_value

    def read(self, path_data=None):

        """
        Load previously built MS/MS spectral index from disk.        
        
        Parameters
        ----------
        path_data : str or Path, optional
            Path to the directory containing the index files.  
            If ``None``, the method uses ``self.path_data``.
            Default is ``None``.

        Returns
        -------
        bool
            ``True`` if the index was successfully read,  
            ``False`` if loading failed for any reason.

        Notes
        -----

        - Any exception during loading returns ``False`` and leaves previous state unchanged.

        """
        

        try:
            if path_data is None:
                path_data = self.path_data

            path_data = Path(path_data)

            with open(path_data / "information_dynamic.json", "r") as f:
                information = json.load(f)

            self.max_ms2_tolerance_in_da = information["max_ms2_tolerance_in_da"]
            self.mass_per_block = information["mass_per_block"]
            self.total_peaks_num = information["total_peaks_num"]
            self.total_spectra_num = information["total_spectra_num"]
            self.extend_fold = information["extend_fold"]
            self.max_indexed_mz = information["max_indexed_mz"]

            is_nl_indexed = information["neutral_loss_indexed"]

            self.index = []
            for name in self.index_names:
                if name == "block_ions_info":
                    self.index.append(np.fromfile(path_data / f"{name}.bin", dtype=self.index_dtypes[name]))

                elif name == "block_nl_info":

                    if is_nl_indexed:
                        self.index.append(np.fromfile(path_data / f"{name}.bin", dtype=self.index_dtypes[name]))
                    else:
                        self.index.append(None)

            return True
        except:
            return False

    def write(self, path_data=None):

        """
        Write the currently built MS/MS spectral index to disk.

        This method serializes all index blocks (product-ion index and optionally the neutral-loss index) along with their associated metadata into the specified directory. 
        The output directory will be created if it does not exist.

        Parameters
        ----------
        path_data : str or Path, optional
            Path to the directory where index files will be saved.  
            If ``None``, the method uses ``self.path_data``.  
            Default is ``None``.

        Returns
        -------
        None
            The method writes files to disk and returns ``None``.

        """
        if path_data is None:
            path_data = self.path_data

        path_data = Path(path_data)
        path_data.mkdir(parents=True, exist_ok=True)
        for i, name in enumerate(self.index_names):
            if self.index[i] is None:
                continue
            self.index[i].tofile(str(path_data / f"{name}.bin"))

        information = {
            "max_ms2_tolerance_in_da": float(self.max_ms2_tolerance_in_da),
            "mass_per_block": float(self.mass_per_block),
            "total_peaks_num": int(self.total_peaks_num),
            "total_spectra_num": int(self.total_spectra_num),
            "extend_fold": float(self.extend_fold),
            "max_indexed_mz": float(self.max_indexed_mz),
            "neutral_loss_indexed": self.index[1] is not None,
        }
        with open(path_data / "information_dynamic.json", "w") as f:
            json.dump(information, f)

        return

    def remove_index(self, path_data=None):

        """
        Remove an existing MS/MS spectral index from disk.

        This method deletes all dynamic index-related files, including block metadata, binary index blocks, and data files. It is typically used after converting an index to a compact format.

        Parameters
        ----------
        path_data : str or Path, optional
            Path to the directory containing the index files to be removed.  
            If ``None``, the method uses ``self.path_data``.  
            Default is ``None``.

        Returns
        -------
        None
            The method removes files from disk and returns ``None``.


        """
        if path_data is None:
            path_data = self.path_data

        path_data = Path(path_data)

        # Remove the block info
        (path_data / "information_dynamic.json").unlink(missing_ok=True)
        for i, name in enumerate(self.index_names):
            file_index = path_data / f"{name}.bin"
            # Delete index
            file_index.unlink(missing_ok=True)

        # Remove the block data
        (path_data / "ions_data.bin").unlink(missing_ok=True)
        (path_data / "nl_data.bin").unlink(missing_ok=True)

        return
    
    def fast_add_new_spectrum_into_index(
        self,
        add_spectrum_list: list,
    ):
        '''
        Add new spectra into the existing index in fast-update mode (without resorting).

        This method appends new spectra to the current on-disk index structure by inserting their fragment ions (and, if available, neutral-loss mass) into existing blocks. 
        Blocks whose reserved capacity is exceeded are moved and expanded, but no global re-sorting of the full index is performed.

        Parameters
        ----------
        add_spectrum_list : list of dict
            A list of spectra to be added. 

        Returns
        -------
        None
            The method updates the on-disk index and in-memory block information in-place.
        
        '''

        # Fast_update mode
        extend_fold = self.extend_fold
        # collect the information of add_spectrum_list
        add_peaks_num = int(np.sum([np.array(spectrum["peaks"]).shape[0] for spectrum in add_spectrum_list]))
        existed_max_spec_idx = self.total_spectra_num - 1
        add_peak_data = self._merge_all_spectra_to_peak_data(add_spectrum_list, add_peaks_num, append=True, existed_max_spec_idx=existed_max_spec_idx)

        [block_ions_info, block_nl_info] = self.index

        ### Open Search ###
        add_peak_data.sort(order="ion_mz")
        search_array = np.arange(0.0, self.max_indexed_mz, self.mass_per_block)
        blocks_num = len(search_array)

        # for open
        with open(self.path_data / "ions_data.bin", "rb+") as f_block_ions_data:
            # build index of add_peak_data for open search
            block_start_loc_array = np.searchsorted(add_peak_data["ion_mz"], search_array, side="left")
            for block_idx in range(blocks_num):
                add_block_start_loc = block_start_loc_array[block_idx]
                if block_idx < blocks_num - 1:
                    add_block_end_loc = block_start_loc_array[block_idx + 1]
                else:
                    add_block_end_loc = len(add_peak_data)

                if add_block_end_loc == add_block_start_loc:
                    continue

                block_start_loc = block_ions_info["start_idx"][block_idx]
                original_part_len = block_ions_info["data_len"][block_idx]
                reserved_len = block_ions_info["reserved_len"][block_idx]

                insert_part_len = add_block_end_loc - add_block_start_loc

                block_data_new_len = original_part_len + insert_part_len
                block_data_insert = np.empty(insert_part_len, dtype=self.dtype_block_data)
                block_data_insert["spec_idx"] = add_peak_data["spec_idx"][add_block_start_loc:add_block_end_loc]
                block_data_insert["mass"] = add_peak_data["ion_mz"][add_block_start_loc:add_block_end_loc]
                block_data_insert["intensity"] = add_peak_data["intensity"][add_block_start_loc:add_block_end_loc]

                block_ions_info["data_len"][block_idx] = block_data_new_len
                block_ions_info["is_sorted"][block_idx] = False
                if block_data_new_len <= reserved_len:
                    insert_location = (block_start_loc + original_part_len) * self.dtype_block_data.itemsize
                    f_block_ions_data.seek(insert_location)
                    block_data_insert.tofile(f_block_ions_data)
                else:
                    # prolong the reserved_len and append to the end and insert the new ion
                    new_reserved_len = int(block_data_new_len * extend_fold)
                    insert_location = f_block_ions_data.seek(0, 2)

                    # Add reserved space for the new block
                    f_block_ions_data.seek(insert_location + new_reserved_len * self.dtype_block_data.itemsize - 1)
                    f_block_ions_data.write(b"1")

                    f_block_ions_data.seek(block_start_loc * self.dtype_block_data.itemsize)
                    original_data = f_block_ions_data.read(original_part_len * self.dtype_block_data.itemsize)

                    f_block_ions_data.seek(insert_location)
                    f_block_ions_data.write(original_data)
                    block_data_insert.tofile(f_block_ions_data)

                    # modify block_info
                    block_ions_info["start_idx"][block_idx] = insert_location // self.dtype_block_data.itemsize
                    block_ions_info["reserved_len"][block_idx] = new_reserved_len

        if block_nl_info is not None:
            add_peak_data.sort(order="nl_mass")
            ### Neutral Loss Search ###
            with open(self.path_data / "nl_data.bin", "rb+") as f_block_nl_data:
                # build index of add_peak_data for neutral_loss search
                block_start_loc_array = np.searchsorted(add_peak_data["nl_mass"], search_array, side="left")
                for block_idx in range(blocks_num):
                    add_block_start_loc = block_start_loc_array[block_idx]
                    if block_idx < blocks_num - 1:
                        add_block_end_loc = block_start_loc_array[block_idx + 1]
                    else:
                        add_block_end_loc = len(add_peak_data)

                    if add_block_end_loc == add_block_start_loc:
                        continue

                    block_start_loc = block_nl_info["start_idx"][block_idx]
                    original_part_len = block_nl_info["data_len"][block_idx]
                    reserved_len = block_nl_info["reserved_len"][block_idx]

                    insert_part_len = add_block_end_loc - add_block_start_loc

                    block_data_insert = np.empty(insert_part_len, dtype=self.dtype_block_data_nl)
                    block_data_insert["spec_idx"] = add_peak_data["spec_idx"][add_block_start_loc:add_block_end_loc]
                    block_data_insert["mass"] = add_peak_data["nl_mass"][add_block_start_loc:add_block_end_loc]
                    block_data_insert["intensity"] = add_peak_data["intensity"][add_block_start_loc:add_block_end_loc]
                    block_data_insert["fragment_mz"] = add_peak_data["ion_mz"][add_block_start_loc:add_block_end_loc]

                    block_nl_info["data_len"][block_idx] = original_part_len + insert_part_len
                    block_nl_info["is_sorted"][block_idx] = False

                    block_data_new_len = original_part_len + insert_part_len
                    if block_data_new_len <= reserved_len:
                        insert_location = (block_start_loc + original_part_len) * self.dtype_block_data_nl.itemsize
                        f_block_nl_data.seek(insert_location)
                        block_data_insert.tofile(f_block_nl_data)
                    else:
                        # prolong the reserved_len and append to the end and insert the new ion
                        new_reserved_len = int(block_data_new_len * extend_fold)
                        insert_location = f_block_nl_data.seek(0, 2)

                        # Add reserved space for the new block
                        f_block_nl_data.seek(insert_location + new_reserved_len * self.dtype_block_data_nl.itemsize - 1)
                        f_block_nl_data.write(b"1")

                        f_block_nl_data.seek(block_start_loc * self.dtype_block_data_nl.itemsize)
                        original_data = f_block_nl_data.read(original_part_len * self.dtype_block_data_nl.itemsize)

                        f_block_nl_data.seek(insert_location)
                        f_block_nl_data.write(original_data)
                        block_data_insert.tofile(f_block_nl_data)

                        # modify block_info
                        block_nl_info["start_idx"][block_idx] = insert_location // self.dtype_block_data_nl.itemsize
                        block_nl_info["reserved_len"][block_idx] = new_reserved_len

        self.index = [block_ions_info, block_nl_info]

        return

    def convert_to_fast_search(
        self,
    ):

        """
        Convert the current index into fast search mode.

        In fast search mode, all peaks stored in the index are sorted by their corresponding mass key:

        - Product-ion index is sorted by fragment m/z.
        - Neutral-loss index (if present) is sorted by neutral-loss mass.


        Parameters
        ----------
        None

        Returns
        -------
        None
            The method updates the on-disk index and in-memory block metadata in-place.

        """
        self.read()
        [block_ions_info, block_nl_info] = self.index
        # Sort open index
        block_ions_data_path = self.path_data / "ions_data.bin"
        original_block_ions_data = np.memmap(block_ions_data_path, dtype=self.dtype_block_data, mode="r+")
        for idx, block_info in enumerate(block_ions_info):

            block_is_sorted = block_info["is_sorted"]

            if block_is_sorted == True:
                continue
            elif block_is_sorted == False:

                block_start_loc = block_info["start_idx"]
                original_part_len = block_info["data_len"]

                block_data = original_block_ions_data[block_start_loc : (block_start_loc + original_part_len)]
                block_data.sort(order="mass")

                block_ions_info[idx]["is_sorted"] = True
        original_block_ions_data.flush()

        if block_nl_info is not None:
            # Sort nl index
            block_nl_data_path = self.path_data / "nl_data.bin"
            original_block_nl_data = np.memmap(block_nl_data_path, dtype=self.dtype_block_data_nl, mode="r+")
            for idx, block_info in enumerate(block_nl_info):

                block_is_sorted = block_info["is_sorted"]

                if block_is_sorted == True:
                    continue
                elif block_is_sorted == False:

                    block_start_loc = block_info["start_idx"]
                    original_part_len = block_info["data_len"]

                    block_data = original_block_nl_data[block_start_loc : (block_start_loc + original_part_len)]
                    block_data.sort(order="mass")

                    block_nl_info[idx]["is_sorted"] = True
            original_block_nl_data.flush()

        self.index = [block_ions_info, block_nl_info]

        return

    def add_new_spectrum_into_index(
        self,
        add_spectrum_list: list,
    ):

        """
        Add new spectra into the index with full sorting.

        This method maintains sorted order within each block by merging existing block data with the new spectra and re-sorting the combined block.

        Parameters
        ----------
        add_spectrum_list : list of dict
            A list of spectra to be added. 

        Returns
        -------
        None
            The method updates the on-disk index and in-memory block metadata in-place.


        """

        extend_fold = self.extend_fold
        # collect the information of add_spectrum_list
        add_peaks_num = int(np.sum([np.array(spectrum["peaks"]).shape[0] for spectrum in add_spectrum_list]))
        existed_max_spec_idx = self.total_spectra_num - 1
        add_peak_data = self._merge_all_spectra_to_peak_data(add_spectrum_list, add_peaks_num, append=True, existed_max_spec_idx=existed_max_spec_idx)

        [
            block_ions_info,
            block_nl_info,
        ] = self.index

        ### Open Search ###
        add_peak_data.sort(order="ion_mz")
        search_array = np.arange(0.0, self.max_indexed_mz, self.mass_per_block)
        blocks_num = len(search_array)

        # for open
        orginal_block_data_open = np.memmap(self.path_data / "ions_data.bin", dtype=self.dtype_block_data, mode="r")
        with open(self.path_data / "ions_data.bin", "rb+") as f_block_ions_data:
            # build index of add_peak_data for open search
            block_start_loc_array = np.searchsorted(add_peak_data["ion_mz"], search_array, side="left")
            for block_idx in range(blocks_num):
                add_block_start_loc = block_start_loc_array[block_idx]
                if block_idx < blocks_num - 1:
                    add_block_end_loc = block_start_loc_array[block_idx + 1]
                else:
                    add_block_end_loc = len(add_peak_data)

                if add_block_end_loc == add_block_start_loc:
                    continue

                block_start_loc = block_ions_info["start_idx"][block_idx]
                original_part_len = block_ions_info["data_len"][block_idx]
                reserved_len = block_ions_info["reserved_len"][block_idx]

                insert_part_len = add_block_end_loc - add_block_start_loc

                block_data_new_len = original_part_len + insert_part_len
                block_data_new = np.empty(original_part_len + insert_part_len, dtype=self.dtype_block_data)
                block_data_new[:original_part_len] = orginal_block_data_open[block_start_loc : (block_start_loc + original_part_len)]
                block_data_new["spec_idx"][original_part_len:] = add_peak_data["spec_idx"][add_block_start_loc:add_block_end_loc]
                block_data_new["mass"][original_part_len:] = add_peak_data["ion_mz"][add_block_start_loc:add_block_end_loc]
                block_data_new["intensity"][original_part_len:] = add_peak_data["intensity"][add_block_start_loc:add_block_end_loc]
                block_data_new.sort(order="mass")
                # if block_data_new_len <= reserved_len, insert directly
                if block_data_new_len <= reserved_len:
                    # insert the new ion
                    insert_location = block_start_loc * self.dtype_block_data.itemsize
                    f_block_ions_data.seek(insert_location)
                    block_data_new.tofile(f_block_ions_data)

                    # modify data_len in block_info
                    block_ions_info["data_len"][block_idx] = original_part_len + insert_part_len

                else:
                    # prolong the reserved_len and append to the end and insert the new ion
                    new_reserved_len = int(block_data_new_len * extend_fold)
                    f_block_ions_data.seek(0, 2)
                    new_start_idx = f_block_ions_data.tell() // self.dtype_block_data.itemsize
                    block_data_new.tofile(f_block_ions_data)
                    # Add reserved space for the new block
                    np.empty(new_reserved_len - block_data_new_len, dtype=self.dtype_block_data).tofile(f_block_ions_data)

                    # modify block_info
                    block_ions_info["start_idx"][block_idx] = new_start_idx
                    block_ions_info["data_len"][block_idx] = block_data_new_len
                    block_ions_info["reserved_len"][block_idx] = new_reserved_len
                block_ions_info["is_sorted"][block_idx] = True

        if block_nl_info is not None:
            add_peak_data.sort(order="nl_mass")
            ### Neutral Loss Search ###
            orginal_block_data_nl = np.memmap(self.path_data / "nl_data.bin", dtype=self.dtype_block_data_nl, mode="r")
            with open(self.path_data / "nl_data.bin", "rb+") as f_block_nl_data:
                # build index of add_peak_data for neutral_loss search
                block_start_loc_array = np.searchsorted(add_peak_data["nl_mass"], search_array, side="left")
                for block_idx in range(blocks_num):
                    add_block_start_loc = block_start_loc_array[block_idx]
                    if block_idx < blocks_num - 1:
                        add_block_end_loc = block_start_loc_array[block_idx + 1]
                    else:
                        add_block_end_loc = len(add_peak_data)

                    if add_block_end_loc == add_block_start_loc:
                        continue

                    block_start_loc = block_nl_info["start_idx"][block_idx]
                    original_part_len = block_nl_info["data_len"][block_idx]
                    reserved_len = block_nl_info["reserved_len"][block_idx]

                    insert_part_len = add_block_end_loc - add_block_start_loc

                    block_data_new_len = original_part_len + insert_part_len
                    block_data_new = np.empty(original_part_len + insert_part_len, dtype=self.dtype_block_data_nl)
                    block_data_new[:original_part_len] = orginal_block_data_nl[block_start_loc : (block_start_loc + original_part_len)]
                    block_data_new["spec_idx"][original_part_len:] = add_peak_data["spec_idx"][add_block_start_loc:add_block_end_loc]
                    block_data_new["mass"][original_part_len:] = add_peak_data["nl_mass"][add_block_start_loc:add_block_end_loc]
                    block_data_new["intensity"][original_part_len:] = add_peak_data["intensity"][add_block_start_loc:add_block_end_loc]
                    block_data_new["fragment_mz"][original_part_len:] = add_peak_data["ion_mz"][add_block_start_loc:add_block_end_loc]
                    block_data_new.sort(order="mass")
                    # if block_data_new_len <= reserved_len, insert directly
                    if block_data_new_len <= reserved_len:
                        # insert the new ion
                        insert_location = block_start_loc * self.dtype_block_data_nl.itemsize
                        f_block_nl_data.seek(insert_location)
                        block_data_new.tofile(f_block_nl_data)

                        # modify data_len in block_info
                        block_nl_info["data_len"][block_idx] = original_part_len + insert_part_len

                    else:
                        # prolong the reserved_len and append to the end and insert the new ion
                        new_reserved_len = int(block_data_new_len * extend_fold)
                        f_block_nl_data.seek(0, 2)
                        new_start_idx = f_block_nl_data.tell() // self.dtype_block_data_nl.itemsize
                        block_data_new.tofile(f_block_nl_data)
                        # Add reserved space for the new block
                        np.empty(new_reserved_len - block_data_new_len, dtype=self.dtype_block_data_nl).tofile(f_block_nl_data)

                        # modify block_info
                        block_nl_info["start_idx"][block_idx] = new_start_idx
                        block_nl_info["data_len"][block_idx] = block_data_new_len
                        block_nl_info["reserved_len"][block_idx] = new_reserved_len
                    block_nl_info["is_sorted"][block_idx] = True

        self.index = [block_ions_info, block_nl_info]

        return

    def get_topn_spec_idx_and_similarity(
        self,
        similarity_array,
        topn=None,
        min_similarity=0.1,
    ):
        '''
        Get the indices and similarity scores of the top-N most similar items.

        This function sorts the similarity array in descending order, selects the top-N indices, and filters out those below the provided minimum similarity threshold.
        
        Parameters
        ----------
        similarity_array : numpy.ndarray
            Array of similarity scores.
        topn : int, optional
            Number of top similarity scores to return. If ``None``, all entries are considered.
        min_similarity : float, optional
            Minimum similarity threshold. Scores below this value are excluded.
            Default is ``0.1``.

        Returns
        -------
        tuple[list[int], list[float]]
            A tuple containing:

            - **result_idx** - List of indices corresponding to selected scores.
            - **result** - List of similarity values for the selected indices.
        '''
        if topn == None:
            topn = len(similarity_array)

        if min_similarity == None:
            min_similarity = 0.0

        # Get topn indices
        topn_indices = np.argsort(similarity_array)[::-1][:topn]

        result = []
        result_idx = []
        for index in topn_indices:
            if similarity_array[index] < min_similarity:
                break

            selected_similarity = similarity_array[index]
            result.append(selected_similarity)
            result_idx.append(index)

        return result_idx, result

    def _extract_data_for_flash(
        self,
    ):
        # Extract data from existing files
        # self.read()
        [block_ions_info, block_nl_info] = self.index

        dtype_ions_flash = np.dtype(
            [
                ("spec_idx", np.uint32),
                ("ion_mz", np.float32),
                ("intensity", np.float32),
            ],
            align=True,
        )

        dtype_nl_flash = np.dtype(
            [
                ("ion_mz", np.float32),
                ("spec_idx", np.uint32),
                ("nl_mass", np.float32),
                ("intensity", np.float32),
            ],
            align=True,
        )

        # For open
        flash_ions = np.zeros(self.total_peaks_num, dtype=dtype_ions_flash)

        block_data_open = np.memmap(self.path_data / "ions_data.bin", mode="r", dtype=self.dtype_block_data)

        data_loc = 0
        for info in block_ions_info:
            block_data_len = info["data_len"]
            if block_data_len == 0:
                continue
            else:
                block_start_loc = info["start_idx"]
                block_data = block_data_open[block_start_loc : (block_start_loc + block_data_len)]

                flash_ions[data_loc : data_loc + block_data_len]["spec_idx"] = block_data["spec_idx"]
                flash_ions[data_loc : data_loc + block_data_len]["ion_mz"] = block_data["mass"]
                flash_ions[data_loc : data_loc + block_data_len]["intensity"] = block_data["intensity"]

                data_loc += block_data_len

        if block_nl_info is not None:
            # For neutral loss
            flash_nl = np.zeros(self.total_peaks_num, dtype=dtype_nl_flash)

            block_data_nl = np.memmap(self.path_data / "nl_data.bin", mode="r", dtype=self.dtype_block_data_nl)

            data_loc = 0
            for info in block_nl_info:
                block_data_len = info["data_len"]
                if block_data_len == 0:
                    continue
                else:
                    block_start_loc = info["start_idx"]
                    block_data = block_data_nl[block_start_loc : (block_start_loc + block_data_len)]

                    flash_nl[data_loc : data_loc + block_data_len]["spec_idx"] = block_data["spec_idx"]
                    flash_nl[data_loc : data_loc + block_data_len]["ion_mz"] = block_data["fragment_mz"]
                    flash_nl[data_loc : data_loc + block_data_len]["nl_mass"] = block_data["mass"]
                    flash_nl[data_loc : data_loc + block_data_len]["intensity"] = block_data["intensity"]

                    data_loc += block_data_len
        else:
            flash_nl = None

        return flash_ions, flash_nl
