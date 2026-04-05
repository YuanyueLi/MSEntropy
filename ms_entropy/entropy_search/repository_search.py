import os
from pathlib import Path

import msgpack
import numpy as np
from ..spectra.tools import clean_spectrum
from ..file_io.spec_file import read_one_spectrum

from .repository_search_without_charge import RepositorySearchWithoutCharge

FILE_METADATA_INFO = np.dtype([("file_id", np.uint64), ("start_loc", np.uint64), ("len", np.uint32)])


class RepositorySearch:
    def __init__(self, path_data: str):

        """
        Initialize :class:`RepositorySearch`.


        Parameters
        ----------
        path_data : str
            Path to the directory containing repository indexes and associated index files.
        """
        self.path_data = Path(path_data)
        self.search_engine: dict[int, RepositorySearchWithoutCharge] = {}
        self._file_metadata_info = None
        self._metadata_index_array = None

    def _get_file_metadata_info(self):
        if self._file_metadata_info is None:
            path_metadata = self.path_data / "metadata"
            os.makedirs(path_metadata, exist_ok=True)
            self._file_metadata_info = open(path_metadata / "metadata_info.bin", "ab+")
        return self._file_metadata_info

    def _save_file_metadata_info(self):
        if self._file_metadata_info is not None:
            self._file_metadata_info.flush()
            os.fsync(self._file_metadata_info.fileno())

    def _get_metadata_index_array(self):
        if self._metadata_index_array is None:
            path_metadata = self.path_data / "metadata"
            try:
                self._metadata_index_array = np.load(path_metadata / "metadata_index.npy")
            except FileNotFoundError:
                os.makedirs(path_metadata, exist_ok=True)
                self._metadata_index_array = np.empty((0,), dtype=FILE_METADATA_INFO)

        return self._metadata_index_array

    def _set_metadata_index_array(self, array: np.ndarray):
        self._metadata_index_array = array

    def _save_metadata_index_array(self):
        path_metadata = self.path_data / "metadata"
        np.save(path_metadata / "metadata_index.npy", self._metadata_index_array)

    def add_ms_file(self, file_name: str, iter_spectra_reader=None):
        """
        Add an MS data file to the repository index.

        This method registers a mass spectrometry data file (e.g. ``mzML``) in the repository metadata index and 
        ingests all MS/MS spectra from the file into the underlying search indices.

        Parameters
        ----------
        file_name : str
            Path to the MS data file (typically an ``mzML`` file).

        iter_spectra_reader : callable, optional
            Iterator function that yields spectrum dictionaries. 
            If ``None`` (default), ``read_one_spectrum`` is used.

        Returns
        -------
        None
        """

        if iter_spectra_reader is None:
            iter_spectra_reader = read_one_spectrum

        # Deal with file metadata
        file_metadata_info = self._get_file_metadata_info()
        metadata_index_array = self._get_metadata_index_array()

        rng = np.random.default_rng()
        file_id = int(rng.integers(low=0, high=np.iinfo(np.int64).max, dtype=np.int64))
        file_metadata = {
            "file_id": file_id,
            "file_name": file_name,
        }
        packed_bytes = msgpack.packb(file_metadata, use_bin_type=True)
        file_metadata_info.seek(0, os.SEEK_END)
        start_loc = file_metadata_info.tell()
        file_metadata_info.write(packed_bytes)
        len_bytes = len(packed_bytes)

        insert_pos = metadata_index_array["file_id"].searchsorted(file_id)
        metadata_index_array = np.insert(metadata_index_array, insert_pos, np.array([(file_id, start_loc, len_bytes)], dtype=FILE_METADATA_INFO))

        self._set_metadata_index_array(metadata_index_array)

        # Deal with spectra
        all_spec = []
        for spec in iter_spectra_reader(file_name):
            if spec["_ms_level"] == 1:
                continue
            spec["scan"] = spec.pop("_scan_number") + file_id
            all_spec.append(spec)

        self._add_new_spectra(spectra_list=all_spec)

    def _add_new_spectra(
        self,
        spectra_list: list,
        insert_mode="fast_update",
        index_for_neutral_loss: bool = True,
        clean_spectra: bool = True,
        precursor_ions_removal_da: float = 1.6,
        noise_threshold=0.01,
        min_ms2_difference_in_da: float = 0.05,
        max_peak_num: int = 0,
    ):

        """
        Preprocess and add spectra to charge-specific search indices.

        This internal method cleans, filters, and partitions spectra by charge state before delegating insertion to the corresponding charge-specific search engines.

        Parameters
        ----------
        spectra_list : list of dict
            List of spectra to be indexed. 

        insert_mode : {"fast_update", "fast_search"}, optional
            Index insertion strategy passed to the underlying charge-specific search engines. 
            Default is ``"fast_update"``.

        index_for_neutral_loss : bool, optional
            Whether to build or update neutral-loss indices in addition to fragment-ion indices. 
            Default is ``True``.

        clean_spectra : bool, optional
            Whether to apply peak cleaning and filtering prior to indexing.
            Default is ``True``.

        precursor_ions_removal_da : float or None, optional
            Precursor ion exclusion window in Dalton. 
            Peaks with m/z values greater than ``precursor_mz - precursor_ions_removal_da`` are removed. 
            If ``None``, no precursor-related filtering is applied.
            Default is ``1.6``.

        noise_threshold : float, optional
            Relative intensity threshold used to remove low-intensity noise peaks during spectrum cleaning. 
            Default is ``0.01``.

        min_ms2_difference_in_da : float, optional
            Minimum allowed m/z difference between adjacent MS/MS peaks after cleaning. 
            Default is ``0.05``.

        max_peak_num : int, optional
            Maximum number of peaks to retain after cleaning. 
            If ``0``, no explicit limit is applied. 
            Default is ``0``.

        Returns
        -------
        None
        """


        spectra_dict = {}
        for spec in spectra_list:
            charge = int(spec["charge"])
            # Clean the peaks
            if abs(charge) == 1:
                if clean_spectra:
                    if precursor_ions_removal_da is not None:
                        max_mz = spec["precursor_mz"] - precursor_ions_removal_da
                    else:
                        max_mz = None

                    spec["peaks"] = clean_spectrum(
                        peaks=spec["peaks"],
                        min_mz=-1,
                        max_mz=max_mz,
                        noise_threshold=noise_threshold,
                        min_ms2_difference_in_da=min_ms2_difference_in_da,
                        max_peak_num=max_peak_num,
                    )

                if len(spec["peaks"]) > 0:
                    if charge not in spectra_dict:
                        spectra_dict[charge] = []
                    spectra_dict[charge].append(spec)

        for charge, specs in spectra_dict.items():
            self._ensure_charge_loaded(charge)
            self.search_engine[charge].add_new_spectra(specs, insert_mode=insert_mode, index_for_neutral_loss=index_for_neutral_loss)

    def build_index(self):
        """
        Build indices for all charge states.

        This method iterates over all loaded charge-specific search engines and triggers index construction for each one. 
        It is typically called after bulk spectrum ingestion.

        Returns
        -------
        None
        """
        for charge, search_engine in self.search_engine.items():
            search_engine.build_index()

    def search_topn_matches(self, method, charge, peaks, precursor_mz=None, topn: int = 1000, output_full_spectrum: bool = False, **kargs):
        
        """
        Search the repository for the top-N matching spectra.

        This method performs a similarity search against the charge-specific search index, and returns ranked match results with optional full-spectrum metadata.

        Parameters
        ----------
        method : str
            Similarity scoring method passed to the underlying search engine.

        charge : int
            Precursor charge state used to select the corresponding charge-specific search index.

        peaks : numpy.ndarray
            Query peak array of shape ``(N, 2)`` containing ``(m/z, intensity)`` pairs.

        precursor_mz : float or None, optional
            Precursor m/z of the query spectrum. Required for certain search modes (e.g. neutral-loss search). 
            Default is ``None``.

        topn : int, optional
            Maximum number of matches to return. 
            Default is ``1000``.

        output_full_spectrum : bool, optional
            If ``True``, return the full indexed spectrum (including peaks and metadata) for each match. 
            If ``False`` (default), only lightweight identifiers and similarity scores are returned.


        Returns
        -------
        list of dict
            List of matching spectra sorted by decreasing similarity. 
        """
        if charge not in self.search_engine:
            charge_dir = self.path_data / f"charge_{charge}"
            if not charge_dir.exists():
                raise ValueError(f"Charge state {charge} not found in the search engine.")
            self.search_engine[charge] = RepositorySearchWithoutCharge(path_data=charge_dir)
            self.search_engine[charge].read()

        result = self.search_engine[charge].search_topn_matches(peaks=peaks, precursor_mz=precursor_mz, method=method, topn=topn, need_metadata=False, **kargs)
        final_results = []

        file_metadata_info = self._get_file_metadata_info()
        metadata_index_array = self._get_metadata_index_array()

        for spec_idx, similarity in result:
            spec_id = self.search_engine[charge].get_spec_id_from_index(spec_idx)
            file_info_idx = metadata_index_array["file_id"].searchsorted(spec_id, side="right") - 1
            file_info_record = metadata_index_array[file_info_idx]
            file_info = msgpack.unpackb(
                os.pread(
                    file_metadata_info.fileno(),
                    file_info_record["len"],
                    file_info_record["start_loc"],
                )
            )
            file_id = file_info["file_id"]

            spec = {
                "file_name": file_info["file_name"],
                "scan": spec_id - file_id,
                "similarity": similarity,
            }

            if output_full_spectrum:
                spec_data = self.search_engine[charge][spec_idx]
                spec_data.update(spec)
                spec = spec_data
            else:
                spec["spec_idx"] = spec_idx

            final_results.append(spec)
        return final_results

    def get_spectrum(self, charge: int, spec_idx: int):

        """
        Retrieve a spectrum by charge state and internal index.

        This method ensures that the charge-specific search engine is loaded and returns the spectrum corresponding to the given internal spectrum index.

        Parameters
        ----------
        charge : int
            Precursor charge state identifying the charge-specific search engine.

        spec_idx : int
            Internal spectrum index within the charge-specific search engine.

        Returns
        -------
        dict
            Spectrum dictionary containing metadata and peak information, as stored in the charge-specific index.
        """
        self._ensure_charge_loaded(charge)
        return self.search_engine[charge][spec_idx]

    def _ensure_charge_loaded(self, charge: int):
        if charge not in self.search_engine:
            charge_dir = self.path_data / f"charge_{charge}"
            self.search_engine[charge] = RepositorySearchWithoutCharge(path_data=charge_dir)

    def write(self):

        """
        Persist repository metadata and search indices to disk.

        This method writes all repository-level metadata structures and charge-specific search indices to their respective search engine instances.

        Returns
        -------
        None
        """
        self._save_file_metadata_info()
        self._save_metadata_index_array()
        for charge, engine in self.search_engine.items():
            engine.write()
