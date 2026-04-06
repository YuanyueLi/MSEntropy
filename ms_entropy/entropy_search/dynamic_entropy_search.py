import pickle
from pathlib import Path

import numpy as np
from .dynamic_entropy_search_core import DynamicEntropySearchCore
from .dynamic_with_flash import DynamicWithFlash
from ..spectra.tools import clean_spectrum


class DynamicEntropySearch:
    def __init__(
        self,
        path_data,
        max_ms2_tolerance_in_da=0.024,
        extend_fold=3,
        mass_per_block: float = 0.05,
        num_per_group: int = 100_000_000,
        cache_list_threshold: int = 1_000_000,
        max_indexed_mz: float = 1500.00005,
        intensity_weight="entropy",  # "entropy" or None
    ):
        """
        Initialize the :class:`DynamicEntropySearch` object.


        Parameters
        ----------
        path_data : str or Path
            Path to the directory where index files are stored.

        max_ms2_tolerance_in_da : float, optional
            Maximum MS/MS tolerance (in Daltons) used during spectrum search.
            Default is ``0.024``.

        extend_fold : int, optional
            Expansion factor for preallocated storage in each m/z block. Determines ``reserved_len = data_len * extend_fold``.  
            Default is ``3``.

        mass_per_block : float, optional
            m/z step size for creating the index blocks.  
            Default is ``0.05 Da``.

        num_per_group : int, optional
            Number of spectra assigned to each group. Default is ``100,000,000``.

        cache_list_threshold : int, optional
            Number of spectra to accumulate in memory before writing them to disk.
            Default is ``1,000,000``.

        max_indexed_mz : float, optional
            Maximum m/z value to index. Ions above this threshold are grouped into a single block. 
            Default is ``1500.00005``.

        intensity_weight : {"entropy", None}, optional
            Determines whether intensities are entropy-weighted.  
            If ``"entropy"``, intensities are weighted accordingly.  
            If ``None``, intensities remain unweighted (equivalent to raw entropy similarity).  
            Default is ``"entropy"``.

        Notes
        -----
        If the index directory already contains ``group_start.pkl`` and ``metadata_start_loc.bin``, they are loaded automatically.
        Otherwise, new metadata structures are initialized.

        The underlying index and search engine is implemented in :class:`DynamicEntropySearchCore`, which is initialized for the most recent group.
        """

        assert cache_list_threshold<=num_per_group, "Cache_list_threshold shouldn't be larger than num_per_group."
        self.path_data = Path(path_data)
        self.num_per_group = num_per_group

        self.max_ms2_tolerance_in_da = max_ms2_tolerance_in_da
        self.extend_fold = extend_fold
        self.mass_per_block = mass_per_block
        self.max_indexed_mz = max_indexed_mz
        self.intensity_weight = intensity_weight

        self.path_data.mkdir(parents=True, exist_ok=True)
        self.cache_list_threshold = cache_list_threshold
        self.cache_list = []
        group_start_path = self.path_data / "group_start.pkl"
        metadata_start_loc_path = self.path_data / "metadata_start_loc.bin"
        if group_start_path.exists() and metadata_start_loc_path.exists():
            self.read()
        else:
            self.group_start = [
                0,
            ]
            self.metadata_start_loc = [
                0,
            ]

        max_group_number = len(self.group_start) - 1
        group_path = self.path_data / f"{max_group_number}"
        self.entropy_search = DynamicEntropySearchCore(
            path_data=group_path,
            max_ms2_tolerance_in_da=self.max_ms2_tolerance_in_da,
            extend_fold=self.extend_fold,
            mass_per_block=self.mass_per_block,
            max_indexed_mz=self.max_indexed_mz,
            intensity_weight=self.intensity_weight,
        )

        self.entropy_search.read()

    def add_new_spectra(
            self, 
            spectra_list: list, 
            insert_mode="fast_update", 
            index_for_neutral_loss: bool = True, 
            convert_to_flash:bool=True, 
            clean=True,
            precursor_ions_removal_da:float=1.6,
            noise_threshold:float=0.01,
            min_ms2_difference_in_da:float=0.05,
            max_peak_num:int=-1
            ):
        
        """
        Add new spectra to the index.

        This function serializes spectrum metadata, appends it to the metadata storage files, updates metadata offsets, and temporarily stores spectra in an in-memory cache. 
        When the cache size exceeds ``cache_list_threshold``, the function triggers incremental index building.

        Parameters
        ----------
        spectra_list : list
            A list of spectra to be added. 
            All spectra must have a single ion mode and should be preprocessed into the correct format.

        insert_mode : {"fast_update", "fast_search"}, optional
            Insertion mode used when updating the index.

            - ``"fast_update"`` (default): Building index structures without resorting.
            - ``"fast_search"``: Building index structures with resorting.

        index_for_neutral_loss : bool, optional
            If ``True`` (default), the index will also maintain entries for neutral-loss ions when building new index blocks.
        
        convert_to_flash : bool, optional
            Whether to convert spectra into a compact format as the :class:`FlashEntropySearch` after the group is full. Default is ``True``.

        clean : bool, optional
            Whether to clean the spectra before adding.  
            Default is ``True``.

        precursor_ions_removal_da : float, optional
            Peaks with m/z greater than ``precursor_mz - precursor_ions_removal_da`` are removed during cleaning. 
            Default is ``1.6`` Da.

        noise_threshold : float, optional
            Relative intensity threshold for noise filtering during cleaning.  
            Peaks with intensity ``< noise_threshold * max(intensity)`` are removed.  
            Default is ``0.01``.

        min_ms2_difference_in_da : float, optional
            Minimum spacing allowed between MS/MS peaks during cleaning.  
            Default is ``0.05`` Da.

        max_peak_num : int or None, optional
            Maximum number of peaks to keep after cleaning.  
            ``None`` keeps all peaks. 
            Default is ``None``.
            
        Notes
        -----

        Spectra are held in ``cache_list`` until the cache reaches the size specified by :attr:`cache_list_threshold`. 
        At that point, :meth:`build_index` is invoked to integrate the cached spectra into the persistent index.

        Returns
        -------
        None
        """

        #Clean
        if clean:
            all_spectra_list=[]

            for spec in spectra_list:

                if 'precursor_mz' not in spec:
                    raise ValueError(f"Spectrum missing 'precursor_mz' field: {spec}")
                if 'peaks' not in spec:
                    raise ValueError(f"Spectrum missing 'peaks' field: {spec}")

                if precursor_ions_removal_da is not None:
                    max_mz = spec['precursor_mz'] - precursor_ions_removal_da

                else:
                    max_mz = -1

                spec['peaks'] = clean_spectrum(
                    peaks=spec['peaks'],
                    min_mz=0,
                    max_mz=max_mz,
                    noise_threshold=noise_threshold,
                    min_ms2_difference_in_da=min_ms2_difference_in_da,
                    max_peak_num=max_peak_num,
                    normalize_intensity=True,
                )

                if len(spec["peaks"]) > 0:
                    all_spectra_list.append(spec)
            
            spectra_list=all_spectra_list

        metadata_array = []
        # Convert metadata of every spectrum to pkl
        for spec in spectra_list:
            metadata_array.append(pickle.dumps(spec))

        metadata_start_loc_path = self.path_data / "metadata_start_loc.bin"
        if metadata_start_loc_path.exists():
            self.metadata_start_loc = np.memmap(metadata_start_loc_path, mode="r", dtype=np.uint64)

        start_loc = [self.metadata_start_loc[-1]]
        if start_loc == [
            0,
        ]:
            metadata_len = np.cumsum(start_loc + [len(metadata) for metadata in metadata_array]).astype(np.uint64)
        else:
            metadata_len = np.cumsum(start_loc + [len(metadata) for metadata in metadata_array]).astype(np.uint64)[1:]

        with open(self.path_data / "metadata_start_loc.bin", "ab") as f:
            metadata_len.tofile(f)

        # Concatenate them
        metadata_all = b"".join(metadata_array)

        # Save the metadata to file
        metadata_path = self.path_data / "metadata.pkl"

        with open(metadata_path, "ab") as f:
            f.write(metadata_all)

        # Move spectra to cache_list
        self.cache_list.extend(spectra_list)

        # Use spectra in cache_list to build index based on the number of spectra in cache_list
        while len(self.cache_list) >= self.cache_list_threshold:

            self.build_index(insert_mode=insert_mode, index_for_neutral_loss=index_for_neutral_loss, convert_to_flash=convert_to_flash)

        return

    def build_index(self, insert_mode="fast_update", index_for_neutral_loss: bool = True, convert_to_flash:bool=True):
        
        """
        Build or update the spectral index.

        This method processes the spectra stored in ``cache_list`` and integrates them into the on-disk index. 
        Depending on the current number of indexed spectra, the method may insert spectra into the existing index, build the index from scratch, or create a new index group.  
        It is also invoked internally by :meth:`add_new_spectra`.

        Parameters
        ----------
        insert_mode : {"fast_update", "fast_search"}, optional
            The index update strategy.

            - ``"fast_update"`` (default): Incrementally appends spectra to the existing index without resorting blocks.
            - ``"fast_search"``: Resorts index blocks for optimized search performance.

        index_for_neutral_loss : bool, optional
            If ``True`` (default), neutral-loss ions are also indexed when building new blocks.

        convert_to_flash : bool, optional
            When a new index group is created, determines whether the existing index is converted to the compact Flash Entropy Search format before writing. 
            Default is ``True``.

        Returns
        -------
        None
        """
        if len(self.cache_list) >= self.cache_list_threshold:
            spectra_to_build = self.cache_list[: self.cache_list_threshold]
            self.cache_list = self.cache_list[self.cache_list_threshold :]

        else:
            if len(self.cache_list) == 0:
                self.entropy_search.write()
                return
            spectra_to_build = self.cache_list
            self.cache_list = []

        # Option 1: if total_spectra_num < self.num_per_group: Insert into an existing group directly
        # Option 2: if total_spectra_num >= self.num_per_group: Create a new group
        # Option 3: if total_spectra_num ==0: build index directly

        # Insert
        if 0 < self.entropy_search.total_spectra_num < self.num_per_group:
            if insert_mode == "fast_update":
                self.entropy_search.fast_add_new_spectrum_into_index(
                    add_spectrum_list=spectra_to_build,
                )

            elif insert_mode == "fast_search":
                # If change to fast_search after fast_update, use `self.convert_to_fast_search()`.
                self.entropy_search.add_new_spectrum_into_index(
                    add_spectrum_list=spectra_to_build,
                )

        # Build directly
        elif self.entropy_search.total_spectra_num == 0:
            self.entropy_search.build_index(all_spectra_list=spectra_to_build, index_for_neutral_loss=index_for_neutral_loss)

        # Create a new class, build de novo
        elif self.entropy_search.total_spectra_num >= self.num_per_group:
            if convert_to_flash:
                # Convert current index to Flash Entropy Search index
                self.entropy_search.convert_to_fast_search()
                self.convert_current_index_to_flash()
            else:
                self.entropy_search.write()

            # Update self.group_start
            self.group_start.append(self.group_start[-1] + self.entropy_search.total_spectra_num)
            max_group_number = len(self.group_start) - 1
            group_path = self.path_data / f"{max_group_number}"
            self.entropy_search = DynamicEntropySearchCore(
                path_data=group_path,
                max_ms2_tolerance_in_da=self.max_ms2_tolerance_in_da,
                extend_fold=self.extend_fold,
                mass_per_block=self.mass_per_block,
                max_indexed_mz=self.max_indexed_mz,
                intensity_weight=self.intensity_weight,
            )

            self.entropy_search.build_index(all_spectra_list=spectra_to_build, index_for_neutral_loss=index_for_neutral_loss)

        # Record precursor_mz
        precursor_mz_array = np.array([spec["precursor_mz"] for spec in spectra_to_build], dtype=np.float32)
        max_group_number = len(self.group_start) - 1
        group_path = self.path_data / f"{max_group_number}"
        with open(group_path / "precursor_mz_array.bin", "ab") as f:
            precursor_mz_array.tofile(f)

        return

    def write(
        self,
    ):
        """
        Write index metadata and group information to disk.

        Notes
        -----
        The following components are written:

        - Internal index data via :meth:`DynamicEntropySearchCore.write`.
        - The list ``group_start`` is serialized to ``group_start.pkl`` to record the cumulative spectrum count at the start of each index group.

        Returns
        -------
        None
        """
        # Write the information of self.entropy_search
        self.entropy_search.write()

        # Write the information of self.group_start
        with open(self.path_data / "group_start.pkl", "wb") as f:
            pickle.dump(self.group_start, f)

        return

    def read(
        self,
    ):
        """
        Load index metadata from disk.

        This method reads previously saved group boundary information and metadata start offsets from the index directory. 

        Notes
        -----
        The following components are loaded:

        - ``group_start`` from ``group_start.pkl``  
        Contains the cumulative spectrum counts marking the start of each index group.

        - ``metadata_start_loc`` from ``metadata_start_loc.bin``  
        A memory-mapped array of byte offsets indicating the starting position of each serialized spectrum metadata entry.

        Returns
        -------
        None
        """

        file_name_group = self.path_data / "group_start.pkl"
        file_name_metadata = self.path_data / "metadata_start_loc.bin"

        with open(file_name_group, "rb") as f:
            self.group_start = pickle.load(f)

        self.metadata_start_loc = np.memmap(file_name_metadata, mode="r", dtype=np.uint64)
        return

    def __getitem__(self, index):

        """
        Retrieve a spectrum's metadata using its global index.

        This method returns the metadata of a spectrum given its global spectrum index within the entire library. 
        The global index corresponds to the cumulative indexing defined by ``group_start`` across all groups.

        Parameters
        ----------
        index : int
            The global spectrum index within the library.  
            This index determines the byte-range in ``metadata.pkl`` using ``metadata_start_loc``.

        Returns
        -------
        dict
            The spectrum metadata corresponding to the requested index.

        """      

        metadata_path = self.path_data / "metadata.pkl"
        self.read()
        with open(metadata_path, "rb") as f:
            f.seek(self.metadata_start_loc[index])
            data = f.read(self.metadata_start_loc[index + 1] - self.metadata_start_loc[index])
            spectrum = pickle.loads(data)

        return spectrum

    def get_metadata(self, group_idx, spec_idx):

        """
        Retrieve the metadata for a spectrum specified by its group and within-group index.

        This method computes the global spectrum index by combining the group-level offset stored in ``group_start`` with the spectrum's position inside that group. 
        It then returns the metadata for the corresponding spectrum.

        Parameters
        ----------
        group_idx : int
            The index of the spectrum group.  
            ``group_start[group_idx]`` gives the global index of the first spectrum in this group.

            For example:
            - Group 0 always begins at global index ``0``.
            - If group 0 contains ``1,000,000`` spectra, then group 1 begins at global index ``1,000,000``; group 2 begins at the cumulative count of groups 0 and 1, etc.

        spec_idx : int
            The index of the spectrum within the specified group.  
            The global spectrum index is computed as:

            ``global_spec_idx = group_start[group_idx] + spec_idx``

        Returns
        -------
        dict
            The metadata dictionary for the requested spectrum.

        Notes
        -----
        This method internally uses :meth:`__getitem__` to retrieve the spectrum metadata once the global index has been determined.

        """
        self.read()
        global_spec_idx = self.group_start[group_idx] + spec_idx
        items = self[global_spec_idx]

        return items

    def identity_search(
        self,
        precursor_mz,
        peaks,
        ms1_tolerance_in_da,
        ms2_tolerance_in_da,
    ):
        """
        Perform an identity search across all indexed spectra.

        This method searches the entire spectral library under ``self.path_data`` for spectra whose precursor m/z values fall within the specified MS1 tolerance of the query precursor. 
        For each matching spectrum, the entropy similarity score is computed using the underlying :class:`DynamicEntropySearchCore` or :class:`DynamicWithFlash` search engine.

        Parameters
        ----------
        precursor_mz : float
            The precursor m/z value of the query spectrum.

        peaks : array-like
            A list or array of peaks from the query spectrum.  
            Peaks must be preprocessed (cleaned, centroided, filtered, etc.) before calling this method.

        ms1_tolerance_in_da : float
            MS1 (precursor m/z) tolerance in Daltons.

        ms2_tolerance_in_da : float
            MS2 (fragment ion) tolerance in Daltons.

        Returns
        -------
        numpy.ndarray
            A 1D array of entropy similarity scores of length ``N``, where ``N`` is the total number of spectra in the library.  
            Spectra whose precursor m/z does not match within tolerance receive a score of ``0.0``.

        """

        all_result = []
        for group in range(len(self.group_start)):
            group_path = self.path_data / f"{group}"
            precursor_mz_file = group_path / "precursor_mz_array.bin"

            if precursor_mz_file.exists():
                precursor_mz_array = np.memmap(precursor_mz_file, mode="r", dtype=np.float32)
                spec_idx = np.where(abs(precursor_mz_array - precursor_mz) <= ms1_tolerance_in_da)[0]
            else:
                raise RuntimeError("Precursor_mz_array not loaded. Call add_new_spectra(...) first. ")
            
            entropy_search=self._assign_entropy_search(group_path=group_path)
                
            cur_result=np.zeros(entropy_search.total_spectra_num, dtype=np.float32)
            cur_result[spec_idx] = entropy_search.search(method="open", peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)[spec_idx]
            all_result.append(cur_result)
            
        identity_result = np.concatenate(all_result)

        return identity_result

    def open_search(
        self,
        peaks,
        ms2_tolerance_in_da,
    ):
        """
        Perform an open search across the entire spectral library.

        This method computes entropy similarity scores between the query spectrum and all spectra stored in the library, using the *open search* strategy of :class:`DynamicEntropySearchCore` or :class:`DynamicWithFlash`. 
        

        Parameters
        ----------
        peaks : array-like
            The peaks of the query spectrum.  
            The spectrum must be preprocessed (e.g., centroided, filtered) before search.

        ms2_tolerance_in_da : float
            Fragment-ion (MS2) tolerance in Daltons.

        Returns
        -------
        numpy.ndarray
            A 1D array of entropy similarity scores with shape ``(N,)``, where ``N`` is the total number of spectra in the library.


        """
        result = []
        # Go through all groups in this library
        for group in range(len(self.group_start)):
            group_path = self.path_data / f"{group}"
            
            entropy_search=self._assign_entropy_search(group_path=group_path)
            
            cur_result = entropy_search.search(method="open", peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)
            result.append(cur_result)

        open_result = np.concatenate(result)
        return open_result

    def neutral_loss_search(self, precursor_mz, peaks, ms2_tolerance_in_da):

        """
        Perform a neutral-loss search across the spectral library.

        This method computes entropy similarity scores based on *neutral-loss* matching, where peaks are compared after subtracting the precursor m/z.

        Parameters
        ----------
        precursor_mz : float
            The precursor m/z value of the query spectrum.

        peaks : array-like
            The cleaned peaks of the query spectrum.  
            Peaks should be preprocessed before calling this method.

        ms2_tolerance_in_da : float
            Fragment-ion (MS2) tolerance in Daltons.

        Returns
        -------
        numpy.ndarray
            A 1D array of entropy similarity scores with shape ``(N,)``, where ``N`` is the total number of spectra in the library.


        """

        result = []
        # Go through all groups in this library
        for group in range(len(self.group_start)):
            group_path = self.path_data / f"{group}"

            entropy_search=self._assign_entropy_search(group_path=group_path)
                
            cur_result = entropy_search.search(method="neutral_loss", precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)
            result.append(cur_result)

        neutral_result = np.concatenate(result)

        return neutral_result

    def hybrid_search(self, precursor_mz, peaks, ms2_tolerance_in_da):

        """
        Perform a hybrid search across the spectral library.

        Hybrid search combines both open search and neutral-loss search to evaluate the similarity between the query spectrum and all spectra in the library. 

        Parameters
        ----------
        precursor_mz : float
            The precursor m/z of the query spectrum.

        peaks : array-like
            The cleaned peaks of the query spectrum.  
            The spectrum should be preprocessed before performing the search.

        ms2_tolerance_in_da : float
            Fragment-ion (MS2) tolerance in Daltons.

        Returns
        -------
        numpy.ndarray
            A 1D array of entropy-based hybrid similarity scores with shape ``(N,)``, where ``N`` is the total number of spectra in the library.

        """
        result = []
        # Go through all groups in this library
        for group in range(len(self.group_start)):
            group_path = self.path_data / f"{group}"

            entropy_search=self._assign_entropy_search(group_path=group_path)

            cur_result = entropy_search.search_hybrid(precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)
            result.append(cur_result)

        hybrid_result = np.concatenate(result)

        return hybrid_result

    def search(
        self,
        precursor_mz,
        peaks,
        ms1_tolerance_in_da=0.01,
        ms2_tolerance_in_da=0.02,
        method="all",
        precursor_ions_removal_da: float = 1.6,
        clean=True,
        noise_threshold=0.01,
        min_ms2_difference_in_da=0.05,
        max_peak_num=None,
    ):
        

        """
        Perform spectral search on a query spectrum.

        This method performs one or more search strategies—including identity, open, neutral-loss, and hybrid search—across the indexed spectral library. 
        The results are returned as a dictionary keyed by search method.

        Parameters
        ----------
        precursor_mz : float
            The precursor m/z value of the query spectrum.

        peaks : array-like
            The MS/MS peaks of the query spectrum, with shape ``(N, 2)`` where each row is ``[mz, intensity]``.  
            Peaks must follow this format and may optionally be cleaned.

        ms1_tolerance_in_da : float, optional
            Tolerance for precursor m/z matching in identity search.
            Default is ``0.01`` Da.

        ms2_tolerance_in_da : float, optional
            Fragment-ion tolerance used by all search modes.
            Default is ``0.02`` Da.

        method : {"identity", "open", "neutral_loss", "hybrid", "all"} or list, optional
            Specifies the search mode(s) to execute.

            - ``"identity"`` - identity search
            - ``"open"`` - open search
            - ``"neutral_loss"`` - neutral-loss search
            - ``"hybrid"`` - combined open + neutral-loss search
            - ``"all"`` (default) - run all four modes


        precursor_ions_removal_da : float, optional
            Peaks with m/z greater than ``precursor_mz - precursor_ions_removal_da`` are removed during spectrum cleaning.  
            Default is ``1.6`` Da.

        clean : bool, optional
            Whether to clean the query spectrum before searching.  
            Default is ``True``.

        noise_threshold : float, optional
            Peaks with intensities below ``noise_threshold * max(intensity)`` are removed during cleaning. 
            Default is ``0.01``.

        min_ms2_difference_in_da : float, optional
            Minimum allowed spacing between MS/MS peaks after cleaning.
            Default is ``0.05`` Da.

        max_peak_num : int or None, optional
            Maximum number of peaks to retain after cleaning.  
            ``None`` (default) keeps all peaks.

        Returns
        -------
        dict
            A dictionary mapping each selected search method to its corresponding entropy similarity score array.  

            Each value is a NumPy array with length equal to the total number of spectra in the library.

        Notes
        -----
        - Cleaning is performed by :func:`clean_spectrum` if ``clean=True``.  
        - When ``method="all"``, all four search strategies are executed.  
        - Only the relevant search modes are invoked based on the ``method`` argument.

        """
        # Assign max_mz
        if precursor_ions_removal_da is not None:
            max_mz = precursor_mz - precursor_ions_removal_da

        else:
            max_mz = -1

        # Clean the query peaks
        if clean:
            peaks = clean_spectrum(
                peaks=peaks,
                min_mz=0,
                max_mz=max_mz,
                noise_threshold=noise_threshold,
                min_ms2_difference_in_da=min_ms2_difference_in_da,
                max_peak_num=max_peak_num,
                normalize_intensity=True,
            )

        # Parse method
        if method == "all":
            method = {"identity", "open", "neutral_loss", "hybrid"}

        elif isinstance(method, str):
            method = {method}

        # Perform search
        result = {}
        if "identity" in method:
            result["identity_search"] = self.identity_search(
                precursor_mz=precursor_mz, peaks=peaks, ms1_tolerance_in_da=ms1_tolerance_in_da, ms2_tolerance_in_da=ms2_tolerance_in_da
            )
        
        if "open" in method:
            result["open_search"] = self.open_search(peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)
        
        if "neutral_loss" in method:
            result["neutral_loss_search"] = self.neutral_loss_search(
                precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da
            )
        
        if "hybrid" in method:
            result["hybrid_search"] = self.hybrid_search(
                precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da
            )

        return result

    def search_topn_matches(
        self,
        peaks,
        precursor_mz=None,
        ms1_tolerance_in_da=0.01,
        ms2_tolerance_in_da=0.02,
        method="open",
        clean=True,
        precursor_ions_removal_da: float = 1.6,
        noise_threshold=0.01,
        min_ms2_difference_in_da=0.05,
        max_peak_num=None,
        topn: int = 3,
        need_metadata: bool = True,
    ):
        
        """
        Search the spectral library and return the top-N most similar spectra.

        This method performs one selected search mode (identity, open, neutral-loss, or hybrid) and returns the best matching spectra according to the similarity scores. 
        Optionally, the query spectrum can be cleaned prior to searching, and metadata for the matched spectra can be returned.

        Parameters
        ----------
        peaks : array-like
            The MS/MS peaks of the query spectrum, with shape ``(N, 2)`` in the format ``[[mz1, intensity1], [mz2, intensity2], ...]``.

        precursor_mz : float, optional
            The precursor m/z of the query spectrum.  
            Required for ``"identity"``, ``"neutral_loss"``, and ``"hybrid"`` modes.

        ms1_tolerance_in_da : float, optional
            MS1 tolerance in Daltons used in identity filtering.  
            Default is ``0.01`` Da.

        ms2_tolerance_in_da : float, optional
            MS2 fragment tolerance in Daltons for similarity computation.  
            Default is ``0.02`` Da.

        method : {"identity", "open", "neutral_loss", "hybrid"}, optional
            The search mode to perform.  
            Default is ``"open"``.

        clean : bool, optional
            Whether to clean the query spectrum before searching.  
            Default is ``True``.

        precursor_ions_removal_da : float, optional
            Peaks with m/z greater than ``precursor_mz - precursor_ions_removal_da`` are removed during cleaning. 
            Default is ``1.6`` Da.

        noise_threshold : float, optional
            Relative intensity threshold for noise filtering during cleaning.  
            Peaks with intensity ``< noise_threshold * max(intensity)`` are removed.  
            Default is ``0.01``.

        min_ms2_difference_in_da : float, optional
            Minimum spacing allowed between MS/MS peaks during cleaning.  
            Default is ``0.05`` Da.

        max_peak_num : int or None, optional
            Maximum number of peaks to keep after cleaning.  
            ``None`` keeps all peaks. 
            Default is ``None``.

        topn : int, optional
            Number of top-matching spectra to return.  
            If ``None``, all spectra are returned.  
            Default is ``3``.

        need_metadata : bool, optional
            If ``True`` (default), return the metadata dictionary for each matched spectrum.  
            If ``False``, return `(global_index, similarity)` tuples instead.
            Default is ``True``.

        Returns
        -------
        list or list of tuples
            If ``need_metadata=True``:  
                A list of metadata dictionaries, each containing the search result and the similarity score under the key.

            If ``need_metadata=False``:  
                A list of tuples ``(global_spec_idx, similarity_score)`` with a one-to-one correspondence between indices and scores.


        """

        all_topn_result = []
        all_topn_result_idx = []

        if method == "identity" or method=="neutral_loss" or method=="hybrid":
            assert precursor_mz is not None, f"Precursor_mz is necessary for {method} search. This parameter should not be None."
        
        # Assign max_mz
        if precursor_ions_removal_da is not None and precursor_mz is not None:
            max_mz = precursor_mz - precursor_ions_removal_da

        else:
            max_mz = -1

        if clean:
            peaks = clean_spectrum(
                peaks=peaks,
                min_mz=0,
                max_mz=max_mz,
                noise_threshold=noise_threshold,
                min_ms2_difference_in_da=min_ms2_difference_in_da,
                max_peak_num=max_peak_num,
                normalize_intensity=True,
            )
            
        for group in range(len(self.group_start)):
            group_path = self.path_data / f"{group}"

            # Perform search
            if method == "identity":

                precursor_mz_file = group_path / "precursor_mz_array.bin"
                if precursor_mz_file.exists():
                    precursor_mz_array = np.memmap(precursor_mz_file, mode="r", dtype=np.float32)
                    spec_idx = np.where(abs(precursor_mz_array - precursor_mz) <= ms1_tolerance_in_da)[0]
                    
                    entropy_search=self._assign_entropy_search(group_path=group_path)

                    result = np.zeros(entropy_search.total_spectra_num, dtype=np.float32)
                    result[spec_idx] = entropy_search.search(method="open", peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)[spec_idx]

                else:
                    raise RuntimeError("Precursor_mz_array not loaded. Call add_new_spectra(...) first. ")

            elif method == "open":

                entropy_search=self._assign_entropy_search(group_path=group_path)

                result = entropy_search.search(method="open", peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)

            elif method == "neutral_loss":
                
                entropy_search=self._assign_entropy_search(group_path=group_path)

                result = entropy_search.search(method="neutral_loss", precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)

            elif method == "hybrid":

                entropy_search=self._assign_entropy_search(group_path=group_path)
                                
                result = entropy_search.search_hybrid(precursor_mz=precursor_mz, peaks=peaks, ms2_tolerance_in_da=ms2_tolerance_in_da)
            
            else:
                raise ValueError(f"Unknown method: {method}. Use 'identity', 'open', 'neutral_loss' or 'hybrid'.")

            # Collect the results
            if topn is None:
                topn=len(result)
            else:
                topn=int(topn)
            topn = min(topn, len(result))
            topn_result_idx, topn_result = entropy_search.get_topn_spec_idx_and_similarity(similarity_array=result, topn=topn)

            # Get the global spec_idx
            topn_result_idx = [self.group_start[group] + idx for idx in topn_result_idx]

            # Combine the results
            all_topn_result.append(topn_result)
            all_topn_result_idx.append(topn_result_idx)

        # Sort and collect topn results
        all_result = np.concatenate(all_topn_result)
        all_result_idx = np.concatenate(all_topn_result_idx)

        collected_idx = np.argsort(all_result)[::-1][:topn]
        collected_result = all_result[collected_idx]
        collected_result_idx = all_result_idx[collected_idx].astype(np.uint64)

        output_result = []
        # Fetch metadata
        if need_metadata == True:
            for i, spec_idx in enumerate(collected_result_idx):
                spec_metadata = self[spec_idx]
                spec_metadata[f"{method}_search_entropy_similarity"] = collected_result[i]
                output_result.append(spec_metadata)

            return output_result

        else:
            # If need_metadata is False, return will be (global_spec_idx, final_result). There will be a consistent one-to-one match between these two elements.
            return list(zip(collected_result_idx, collected_result))

    def _assign_entropy_search(
            self,
            group_path:Path
    ):
        
        """
        Assign and initialize the appropriate entropy search engine for a given group.

        This method inspects the index directory for the specified group and determines whether the group uses the dynamic index format (``information_dynamic.json``) or 
        the Flash Entropy Search format (``information.json``). It then loads and returns an initialized search engine instance corresponding to the detected format.

        Parameters
        ----------
        group_path : Path
            The path to the directory containing the index files for a specific group of spectra.

        Returns
        -------
        object
            An initialized entropy search object.  


        """
        
        if (group_path/"information_dynamic.json").exists():
            entropy_search=self.entropy_search
            entropy_search.path_data=group_path
        elif (group_path/"information.json").exists():
            entropy_search = DynamicWithFlash(
                path_data=group_path, max_ms2_tolerance_in_da=self.max_ms2_tolerance_in_da, intensity_weight=self.intensity_weight
            )
        else:
            raise FileNotFoundError("Neither information.json nor information_dynamic.json exists. Failed to assign entropy search.")
        
        entropy_search.read()

        return entropy_search
        
    def convert_to_fast_search(
        self,
    ):
        """
        Convert all dynamic index groups to the Flash Entropy Search format.

        This method iterates through all index groups in the library and converts any group stored in dynamic indexing format (identified by the presence of ``information_dynamic.json``) 
        into the compact and search-optimized *Flash Entropy Search* format.

        This operation should be performed **after** calling :meth:`write`, ensuring that all index data is fully written before conversion.


        Returns
        -------
        None
        """

        # If this function is to be used, use it after `write()`.
        for group in range(len(self.group_start)):
            group_path = self.path_data / f"{group}"

            if (group_path/"information_dynamic.json").exists():
                self.entropy_search.path_data = group_path
                self.entropy_search.read()
                self.entropy_search.convert_to_fast_search()
                self.entropy_search.write()

        return

    def convert_current_index_to_flash(self):
        """
        Convert the current dynamic index into the Flash Entropy Search format.

        This method extracts the peak data from the currently active dynamic index and rebuilds it as a Flash-format index using :class:`DynamicWithFlash`.
        After conversion, the original dynamic index files are removed and replaced with the Flash-formatted index.

        This method is intended to be used **after** calling :meth:`convert_to_fast_search`, which ensures that all groups except the active one have already been converted.


        Returns
        -------
        None

        """
        # Use this function after using `convert_to_fast_search()`
        flash_ions, flash_nl = self.entropy_search._extract_data_for_flash()
        dynamic_with_flash = DynamicWithFlash(
            path_data=self.entropy_search.path_data, max_ms2_tolerance_in_da=self.max_ms2_tolerance_in_da, intensity_weight=self.intensity_weight
        )

        dynamic_with_flash.total_peaks_num = self.entropy_search.total_peaks_num
        dynamic_with_flash.total_spectra_num = self.entropy_search.total_spectra_num
        dynamic_with_flash.max_ms2_tolerance_in_da = self.entropy_search.max_ms2_tolerance_in_da

        dynamic_with_flash.build_index(peak_data=flash_ions, max_indexed_mz=self.max_indexed_mz, peak_data_nl=flash_nl)
        self.entropy_search.remove_index()
        dynamic_with_flash.write()
        return
