#!/usr/bin/env python3
import numpy as np
import pickle
from pathlib import Path
from .flash_entropy_search_core import FlashEntropySearchCore
from .flash_entropy_search_core_low_memory import FlashEntropySearchCoreLowMemory
from .flash_entropy_search_core_medium_memory import FlashEntropySearchCoreMediumMemory
from ..spectra import clean_spectrum


class FlashEntropySearch:
    def __init__(
        self,
        max_ms2_tolerance_in_da=0.024,
        mz_index_step=0.0001,
        low_memory=False,
        path_data=None,
        intensity_weight="entropy",
    ):
        """
        Initialize the EntropySearch class.
        :param max_ms2_tolerance_in_da:  The maximum MS2 tolerance in Da.
        :param mz_index_step:   The step size for the m/z index.
        :param low_memory:  The memory usage mode, can be 0, 1, or 2. 0 means normal mode, 1 means low memory mode, and 2 means medium memory mode.
        :param path_data:   The path to save the index data.
        :param intensity_weight:    The weight for the intensity in the entropy calculation, can be "entropy" or None. Default is "entropy".
            - None: The intensity will not be weighted, then the unweighted similarity will be calculated.
            - "entropy": The intensity will be weighted by the entropy, then the entropy similarity will be calculated.
        """
        self.precursor_mz_array = np.zeros(0, dtype=np.float32)
        self.low_memory = low_memory
        if low_memory == 1:
            self.entropy_search = FlashEntropySearchCoreLowMemory(
                path_data=path_data, max_ms2_tolerance_in_da=max_ms2_tolerance_in_da, mz_index_step=mz_index_step, intensity_weight=intensity_weight
            )
        elif low_memory == 2:
            self.entropy_search = FlashEntropySearchCoreMediumMemory(
                path_data=path_data, max_ms2_tolerance_in_da=max_ms2_tolerance_in_da, mz_index_step=mz_index_step, intensity_weight=intensity_weight
            )
        else:
            self.entropy_search = FlashEntropySearchCore(
                path_data=path_data, max_ms2_tolerance_in_da=max_ms2_tolerance_in_da, mz_index_step=mz_index_step, intensity_weight=intensity_weight
            )

    def identity_search(self, precursor_mz, peaks, ms1_tolerance_in_da, ms2_tolerance_in_da, target="cpu", output_matched_peak_number=False, **kwargs):
        """
        Run the identity search, the query spectrum should be preprocessed by `clean_spectrum()` function before calling this function.

        For super large spectral library, directly identity search is not recommended. To do the identity search on super large spectral library,
        divide the spectral library into several parts, build the index for each part, and then do the identity search on each part will be much faster.

        :param precursor_mz:    The precursor m/z of the query spectrum.
        :param peaks:           The peaks of the query spectrum, should be the output of `clean_spectrum()` function.
        :param ms1_tolerance_in_da:  The MS1 tolerance in Da.
        :param ms2_tolerance_in_da:  The MS2 tolerance in Da.
        :param target:  The target device for the search, can be "cpu" or "gpu".
        :param output_matched_peak_number:  If True, the number of matched peaks will be returned with the entropy similarity score.

        :return:    The entropy similarity score for each spectrum in the library, a numpy array with shape (N,), N is the number of spectra in the library.
                    If `output_matched_peak_number` is True, the number of matched peaks will be returned with the entropy similarity score, i.e. the return
                    will be a tuple of two numpy arrays, the first one is the entropy similarity score, and the second one is the number of matched peaks.
        """
        precursor_mz_min = precursor_mz - ms1_tolerance_in_da
        precursor_mz_max = precursor_mz + ms1_tolerance_in_da
        spectra_idx_min = np.searchsorted(self.precursor_mz_array, precursor_mz_min, side="left")
        spectra_idx_max = np.searchsorted(self.precursor_mz_array, precursor_mz_max, side="right")
        if spectra_idx_min >= spectra_idx_max:
            if output_matched_peak_number:
                return np.zeros(self.entropy_search.total_spectra_num, dtype=np.float32), np.zeros(self.entropy_search.total_spectra_num, dtype=np.uint16)
            else:
                return np.zeros(self.entropy_search.total_spectra_num, dtype=np.float32)
        else:
            return self.entropy_search.search(
                method="open",
                target=target,
                peaks=peaks,
                ms2_tolerance_in_da=ms2_tolerance_in_da,
                search_type=1,
                search_spectra_idx_min=spectra_idx_min,
                search_spectra_idx_max=spectra_idx_max,
                output_matched_peak_number=output_matched_peak_number,
            )

    def open_search(self, peaks, ms2_tolerance_in_da, target="cpu", output_matched_peak_number=False, **kwargs):
        """
        Run the open search, the query spectrum should be preprocessed by `clean_spectrum()` function before calling this function.

        :param peaks:           The peaks of the query spectrum, should be the output of `clean_spectrum()` function.
        :param ms2_tolerance_in_da:  The MS2 tolerance in Da.
        :param target:  The target device for the search, can be "cpu" or "gpu".
        :param output_matched_peak_number:  If True, the number of matched peaks will be returned with the entropy similarity score.

        :return:    The entropy similarity score for each spectrum in the library, a numpy array with shape (N,), N is the number of spectra in the library.
                    If `output_matched_peak_number` is True, the number of matched peaks will be returned with the entropy similarity score, i.e. the return
                    will be a tuple of two numpy arrays, the first one is the entropy similarity score, and the second one is the number of matched peaks.
        """
        return self.entropy_search.search(
            method="open",
            target=target,
            peaks=peaks,
            ms2_tolerance_in_da=ms2_tolerance_in_da,
            search_type=0,
            output_matched_peak_number=output_matched_peak_number,
        )

    def neutral_loss_search(self, precursor_mz, peaks, ms2_tolerance_in_da, target="cpu", output_matched_peak_number=False, **kwargs):
        """
        Run the neutral loss search, the query spectrum should be preprocessed by `clean_spectrum()` function before calling this function.

        :param precursor_mz:    The precursor m/z of the query spectrum.
        :param peaks:           The peaks of the query spectrum, should be the output of `clean_spectrum()` function.
        :param ms2_tolerance_in_da:  The MS2 tolerance in Da.
        :param target:  The target device for the search, can be "cpu" or "gpu".
        :param output_matched_peak_number:  If True, the number of matched peaks will be returned with the entropy similarity score.

        :return:    The entropy similarity score for each spectrum in the library, a numpy array with shape (N,), N is the number of spectra in the library.
                    If `output_matched_peak_number` is True, the number of matched peaks will be returned with the entropy similarity score, i.e. the return
                    will be a tuple of two numpy arrays, the first one is the entropy similarity score, and the second one is the number of matched peaks.
        """
        return self.entropy_search.search(
            method="neutral_loss",
            target=target,
            precursor_mz=precursor_mz,
            peaks=peaks,
            ms2_tolerance_in_da=ms2_tolerance_in_da,
            search_type=0,
            output_matched_peak_number=output_matched_peak_number,
        )

    def hybrid_search(self, precursor_mz, peaks, ms2_tolerance_in_da, target="cpu", **kwargs):
        """
        Run the hybrid search, the query spectrum should be preprocessed by `clean_spectrum()` function before calling this function.

        :param precursor_mz:    The precursor m/z of the query spectrum.
        :param peaks:           The peaks of the query spectrum, should be the output of `clean_spectrum()` function.
        :param ms2_tolerance_in_da:  The MS2 tolerance in Da.
        :param target:  The target device for the search, can be "cpu" or "gpu".

        :return:    The entropy similarity score for each spectrum in the library, a numpy array with shape (N,), N is the number of spectra in the library.
        """
        return self.entropy_search.search_hybrid(target=target, precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)

    def clean_spectrum_for_search(
        self, precursor_mz, peaks, precursor_ions_removal_da: float = 1.6, noise_threshold=0.01, min_ms2_difference_in_da: float = 0.05, max_peak_num: int = 0
    ):
        """
        Clean the MS/MS spectrum, need to be called before any search.

        :param precursor_mz:    The precursor m/z of the spectrum.
        :param peaks:           The peaks of the spectrum, should be a list or numpy array with shape (N, 2), N is the number of peaks. The format of the peaks is [[mz1, intensity1], [mz2, intensity2], ...].
        :param precursor_ions_removal_da:   The ions with m/z larger than precursor_mz - precursor_ions_removal_da will be removed.
                                            Default is 1.6. Set to None to not remove any ions.
        :param noise_threshold: The intensity threshold for removing the noise peaks. The peaks with intensity smaller than noise_threshold * max(intensity)
                                will be removed. Default is 0.01.
        :param min_ms2_difference_in_da:    The minimum difference between two peaks in the MS/MS spectrum. Default is 0.05.
        :param max_peak_num:    The maximum number of peaks in the MS/MS spectrum. Default is 0, which means no limit.
        """
        if precursor_ions_removal_da is not None:
            max_mz = precursor_mz - precursor_ions_removal_da
        else:
            max_mz = None
        return clean_spectrum(
            peaks=peaks,
            min_mz=-1,
            max_mz=max_mz,
            noise_threshold=noise_threshold,
            min_ms2_difference_in_da=min_ms2_difference_in_da,
            max_peak_num=max_peak_num,
            normalize_intensity=True,
        )

    def search(
        self,
        precursor_mz,
        peaks,
        ms1_tolerance_in_da=0.01,
        ms2_tolerance_in_da=0.02,
        method="all",  # "identity", "open", "neutral_loss", "hybrid", "all", or list of the above
        target="cpu",
        precursor_ions_removal_da: float = 1.6,
        noise_threshold=0.01,
        min_ms2_difference_in_da: float = 0.05,
        max_peak_num: int = None,
    ):
        """
        Run the Flash entropy search for the query spectrum.

        :param precursor_mz:    The precursor m/z of the query spectrum.
        :param peaks:           The peaks of the query spectrum, should be a list or numpy array with shape (N, 2), N is the number of peaks. The format of the peaks is [[mz1, intensity1], [mz2, intensity2], ...].
        :param ms1_tolerance_in_da:  The MS1 tolerance in Da. Default is 0.01.
        :param ms2_tolerance_in_da:  The MS2 tolerance in Da. Default is 0.02.
        :param method:  The search method, can be "identity", "open", "neutral_loss", "hybrid", "all", or list of the above.
        :param target:  The target device for the search, can be "cpu" or "gpu".
        :param precursor_ions_removal_da:   The ions with m/z larger than precursor_mz - precursor_ions_removal_da will be removed.
                                            Default is 1.6.
        :param noise_threshold: The intensity threshold for removing the noise peaks. The peaks with intensity smaller than noise_threshold * max(intensity)
                                will be removed. Default is 0.01.
        :param min_ms2_difference_in_da:    The minimum difference between two peaks in the MS/MS spectrum. Default is 0.05.
        :param max_peak_num:    The maximum number of peaks in the MS/MS spectrum. Default is None, which means no limit.

        :return:    A dictionary with the search results. The keys are "identity_search", "open_search", "neutral_loss_search", "hybrid_search", and the values are the search results for each method.
        """
        if precursor_ions_removal_da is not None:
            max_mz = precursor_mz - precursor_ions_removal_da
        else:
            max_mz = -1
        peaks = clean_spectrum(
            peaks=peaks,
            min_mz=0,
            max_mz=max_mz,
            noise_threshold=noise_threshold,
            min_ms2_difference_in_da=min_ms2_difference_in_da,
            max_peak_num=max_peak_num,
            normalize_intensity=True,
        )
        if method == "all":
            method = {"identity", "open", "neutral_loss", "hybrid"}
        elif isinstance(method, str):
            method = {method}

        result = {}
        if "identity" in method:
            result["identity_search"] = self.identity_search(
                precursor_mz=precursor_mz, peaks=peaks, ms1_tolerance_in_da=ms1_tolerance_in_da, ms2_tolerance_in_da=ms2_tolerance_in_da, target=target
            )
        if "open" in method:
            result["open_search"] = self.open_search(peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da, target=target)
        if "neutral_loss" in method:
            result["neutral_loss_search"] = self.neutral_loss_search(
                precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da, target=target
            )
        if "hybrid" in method:
            result["hybrid_search"] = self.hybrid_search(precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da, target=target)
        return result

    def build_index(
        self,
        all_spectra_list: list = None,
        max_indexed_mz: float = 1500.00005,
        precursor_ions_removal_da: float = 1.6,
        noise_threshold=0.01,
        min_ms2_difference_in_da: float = 0.05,
        max_peak_num: int = 0,
        clean_spectra: bool = True,
    ):
        """
        Set the library spectra for entropy search.

        The `all_spectra_list` must be a list of dictionaries, with each dictionary containing at least two keys: "precursor_mz" and "peaks".
        The dictionary should be in the format of {"precursor_mz": precursor_mz, "peaks": peaks, ...}, All keys in the dictionary, except "peaks,"
        will be saved as the metadata and can be accessed using the  __getitem__ function (e.g. entropy_search[0] returns the metadata for the
        first spectrum in the library).

            - The precursor_mz is the precursor m/z value of the MS/MS spectrum;

            - The peaks is an numpy array or a nested list of the MS/MS spectrum, looks like [[mz1, intensity1], [mz2, intensity2], ...].

        :param all_spectra_list:    A list of dictionaries in the format of {"precursor_mz": precursor_mz, "peaks": peaks},
                                    the spectra in the list do not need to be sorted by the precursor m/z.
                                    This function will sort the spectra by the precursor m/z and output the sorted spectra list.

        :param max_indexed_mz: The maximum m/z value that will be indexed. Default is 1500.00005.
        :param precursor_ions_removal_da:   The ions with m/z larger than precursor_mz - precursor_ions_removal_da will be removed.
                                            Default is 1.6. Set to None to not remove any ions.
        :param noise_threshold: The intensity threshold for removing the noise peaks. The peaks with intensity smaller than noise_threshold * max(intensity)
                                will be removed. Default is 0.01.
        :param min_ms2_difference_in_da:    The minimum difference between two peaks in the MS/MS spectrum. Default is 0.05.
        :param max_peak_num:    The maximum number of peaks in the MS/MS spectrum. Default is 0, which means no limit.
        :param clean_spectra:   If True, the spectra will be cleaned before indexing. Default is True. If ALL spectra in the library are pre-cleaned with the
                                function `clean_spectrum` or `clean_spectrum_for_search`, set this parameter to False. ALWAYS set this parameter to true if
                                the spectra are not pre-prepossessed with the function `clean_spectrum` or `clean_spectrum_for_search`.

        :return:    If the all_spectra_list is provided, this function will return the sorted spectra list.
        """

        # Sort the spectra by the precursor m/z.
        all_sorted_spectra_list = sorted(all_spectra_list, key=lambda x: x["precursor_mz"])

        # Clean the spectra, and collect the non-empty spectra
        all_spectra_list = []
        all_metadata_list = []
        for spec in all_sorted_spectra_list:
            # Clean the peaks
            if clean_spectra:
                spec["peaks"] = self.clean_spectrum_for_search(
                    peaks=spec["peaks"],
                    precursor_mz=spec["precursor_mz"],
                    precursor_ions_removal_da=precursor_ions_removal_da,
                    noise_threshold=noise_threshold,
                    min_ms2_difference_in_da=min_ms2_difference_in_da,
                    max_peak_num=max_peak_num,
                )
            if len(spec["peaks"]) > 0:
                all_spectra_list.append(spec)
                all_metadata_list.append(pickle.dumps(spec))

        # Extract precursor m/z array
        self.precursor_mz_array = np.array([spec["precursor_mz"] for spec in all_spectra_list], dtype=np.float32)

        # Extract metadata array
        all_metadata_len = np.array([0] + [len(metadata) for metadata in all_metadata_list], dtype=np.uint64)
        self.metadata_loc = np.cumsum(all_metadata_len).astype(np.uint64)
        self.metadata = np.frombuffer(b"".join(all_metadata_list), dtype=np.uint8)

        # Call father class to build the index.
        self.entropy_search.build_index(all_spectra_list, max_indexed_mz)
        return all_spectra_list

    def __getitem__(self, index):
        """
        Get the MS/MS metadate by the index.

        :param index:   The index of the MS/MS spectrum.
        :return:    The MS/MS spectrum in the format of {"precursor_mz": precursor_mz, "peaks": peaks}.
        """
        if self.metadata is not None:
            spectrum = pickle.loads(self.metadata[self.metadata_loc[index] : self.metadata_loc[index + 1]].tobytes())
        else:
            spectrum = {"precursor_mz": self.precursor_mz_array[index]}
        return spectrum

    def get_topn_matches(self, similarity_array, topn=3, min_similarity=0.01):
        """
        Get the topn MS/MS spectra with the highest entropy similarity.

        :param similarity_array:  The entropy similarity of the MS/MS spectra.
        :param topn:    The number of MS/MS spectra to return, if None, all the MS/MS spectra will be returned.
        :param min_similarity:  The minimum similarity of the MS/MS spectra to return, if None, all the MS/MS spectra will be returned.
        :return:    The topn MS/MS spectra with the highest entropy similarity.
        """
        if topn is None:
            topn = len(similarity_array)
        if min_similarity is None:
            min_similarity = 0.0

        # Get the topn indices
        topn_indices = np.argsort(similarity_array)[::-1][:topn]

        result = []
        for index in topn_indices:
            if similarity_array[index] < min_similarity:
                break
            item = self[index]
            item["entropy_similarity"] = similarity_array[index]
            result.append(item)
        return result

    def write(self, path_data=None):
        """
        Write the MS/MS spectral library to a file.

        :param path_data:   The path of the file to write.
        :return:    None
        """
        if path_data is None:
            path_data = self.entropy_search.path_data
        else:
            path_data = Path(path_data)

        path_data = Path(path_data)
        path_data.mkdir(parents=True, exist_ok=True)

        self.precursor_mz_array.tofile(str(path_data / "precursor_mz.npy"))
        self.metadata.tofile(str(path_data / "metadata.npy"))
        self.metadata_loc.tofile(str(path_data / "metadata_loc.npy"))

        self.entropy_search.write(path_data)

    def read(self, path_data=None):
        """
        Read the MS/MS spectral library from a file.

        :param path_data:   The path of the file to read.
        :return:    None
        """
        if path_data is None:
            path_data = self.entropy_search.path_data
        else:
            path_data = Path(path_data)

        if self.low_memory:
            self.precursor_mz_array = np.memmap(path_data / "precursor_mz.npy", dtype=np.float32, mode="r")
            self.metadata = np.memmap(path_data / "metadata.npy", dtype=np.uint8, mode="r")
            self.metadata_loc = np.memmap(path_data / "metadata_loc.npy", dtype=np.uint64, mode="r")
        else:
            self.precursor_mz_array = np.fromfile(str(path_data / "precursor_mz.npy"), dtype=np.float32)
            self.metadata = np.fromfile(str(path_data / "metadata.npy"), dtype=np.uint8)
            self.metadata_loc = np.fromfile(str(path_data / "metadata_loc.npy"), dtype=np.uint64)

        return self.entropy_search.read(path_data)

    def save_memory_for_multiprocessing(self):
        """
        Save the memory for multiprocessing. This function will move the numpy array in the index to shared memory in order to save memory.

        This function is not required when you only use one thread to search the MS/MS spectra.
        When use multiple threads, this function is also not required but highly recommended, as it avoids the memory copy and saves a lot of memory and time.

        :return:    None
        """
        self.entropy_search.save_memory_for_multiprocessing()
