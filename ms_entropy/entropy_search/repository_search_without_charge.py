import numpy as np

from .dynamic_entropy_search import DynamicEntropySearch
from .ms_spectrum import read_spectrum_from_file_stream, write_spectrum_to_file_stream


class RepositorySearchWithoutCharge(DynamicEntropySearch):
    def __init__(self, path_data: str, *args, **kwargs):

        """
        Initialize :class:`RepositorySearchWithoutCharge`.

        Parameters
        ----------
        path_data : str
            Path to the directory containing index files. 
            The directory is created and managed by the parent class.

        """
        super().__init__(path_data, *args, **kwargs)

        self.file_spec_id_array = open(self.path_data / "spec_id_array.bin", "ab")
        self.file_loc_array = open(self.path_data / "loc_array.bin", "ab")

        self.spec_id_array = None
        self.loc_array = None
        self.spec_metadata = None

    def __getitem__(self, index):
        if self.loc_array is None or self.spec_metadata is None:
            self.read()

        self.spec_metadata.seek(self.loc_array[index])
        spec = next(read_spectrum_from_file_stream(self.spec_metadata))
        return spec

    def get_spec_id_from_index(self, index):

        """
        Retrieve a spectrum identifier by its index position.

        This method returns the spectrum ID corresponding to the given index in the spectrum ID array. 
        If the array has not yet been loaded into memory, it is read from disk.

        Parameters
        ----------
        index : int or array-like
            Index of the spectrum to retrieve.

        Returns
        -------
        int or numpy.ndarray
            Spectrum identifier corresponding to the provided index.

        """
        if self.spec_id_array is None:
            self.read()
        return self.spec_id_array[index]

    def add_new_spectra(self, spectra_list: list, insert_mode="fast_update", index_for_neutral_loss: bool = True):

        """
        Add new spectra to index.

        This method appends spectra metadata and binary spectrum data to disk, records spectrum identifiers and byte offsets, 
        and places the spectra into an in-memory cache. 
        When the cache size reaches the configured threshold, the index is built according to the specified insertion mode.

        Parameters
        ----------
        spectra_list : list of dict
            List of spectra to be added. Each spectrum must be a dictionary.

        insert_mode : {"fast_update", "fast_search"}, optional
            Strategy used when building the index:

            -``"fast_update"`` - Building without resorting.

            -``"fast_search"`` - Building with resorting.

            Default is ``"fast_update"``.

        index_for_neutral_loss : bool, optional
            Whether to build or update the neutral-loss index in addition to the fragment-ion index. 

            Default is ``True``.

        Returns
        -------
        None

        """
        
        scan_array = np.array([spec["scan"] for spec in spectra_list], dtype=np.uint64)
        scan_array.tofile(self.file_spec_id_array)

        loc_array = np.empty(len(spectra_list), dtype=np.uint64)
        with open(self.path_data / "ms2_data.bin", "ab") as f_spec_metadata:
            for i, spec in enumerate(spectra_list):
                loc_array[i] = f_spec_metadata.tell()
                write_spectrum_to_file_stream(spec, f_spec_metadata)
        loc_array.tofile(self.file_loc_array)

        # Move spectra to cache_list
        self.cache_list.extend(spectra_list)

        # Use spectra in cache_list to build index based on the number of spectra in cache_list
        while len(self.cache_list) >= self.cache_list_threshold:
            self.build_index(insert_mode=insert_mode, index_for_neutral_loss=index_for_neutral_loss)
        return

    def write(self):
        """
        Flush index state.

        Returns
        -------
        None
        """
        self.file_spec_id_array.close()
        self.file_spec_id_array = open(self.path_data / "spec_id_array.bin", "ab")

        self.file_loc_array.close()
        self.file_loc_array = open(self.path_data / "loc_array.bin", "ab")

        self.entropy_search.write()
        with open(self.path_data / "group_start.bin", "wb") as f:
            np.array(self.group_start, dtype=np.uint64).tofile(f)

    def read(self):

        """
        Load on-disk index structures into memory.

        This method initializes memory-mapped arrays and file handles required for spectrum lookup and random access. 

        Returns
        -------
        None
        """

        self.spec_id_array = np.memmap(self.path_data / "spec_id_array.bin", mode="r", dtype=np.uint64)
        self.loc_array = np.memmap(self.path_data / "loc_array.bin", mode="r", dtype=np.uint64)
        self.spec_metadata = open(self.path_data / "ms2_data.bin", "rb")

        with open(self.path_data / "group_start.bin", "rb") as f:
            self.group_start = np.fromfile(f, dtype=np.uint64).tolist()
