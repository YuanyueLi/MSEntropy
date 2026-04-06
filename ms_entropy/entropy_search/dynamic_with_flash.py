from .flash_entropy_search_core_for_dynamic_indexing import FlashEntropySearchCoreForDynamicIndexing
import numpy as np

class DynamicWithFlash(FlashEntropySearchCoreForDynamicIndexing):
    def build_index(self, peak_data, max_indexed_mz, peak_data_nl=None):

        """
        Build the Flash-format fragment-ion and neutral-loss index.

        This method constructs the index structures required for fast entropy search using the Flash Entropy Search representation. 
        Fragment ions are indexed by their m/z values, and neutral-loss ions (if provided) are indexed by their neutral-loss mass.

        Parameters
        ----------
        peak_data : structured array
            A structured NumPy array containing fragment-ion information.

        max_indexed_mz : float
            Maximum m/z boundary.

        peak_data_nl : structured array, optional
            Structured NumPy array describing neutral-loss peaks:

            If ``None`` (default), no neutral-loss index is constructed.

        Returns
        -------
        tuple
            A tuple ``(fragment_ion_index, neutral_loss_index)``, where:

            - ``fragment_ion_index`` is a tuple containing ``(mz_idx_start, ion_mz, intensity, spec_idx)``

            - ``neutral_loss_index`` is either a tuple containing ``(nl_idx_start, nl_mass, intensity, spec_idx, ion_mz)``, or ``None`` if no neutral-loss data is provided.


        """
        # Record the m/z, intensity, and spectrum index information for product ions.
        all_ions_mz = peak_data["ion_mz"]
        all_ions_intensity = peak_data["intensity"]
        all_ions_spec_idx = peak_data["spec_idx"]

        # Build index for fast access to the ion's m/z.
        all_ions_mz_idx_start = self._generate_index(all_ions_mz, max_indexed_mz)
        fragment_ion_index = (
            all_ions_mz_idx_start,
            all_ions_mz,
            all_ions_intensity,
            all_ions_spec_idx,
        )

        if peak_data_nl is not None:

            # Build the index for fast access to the neutral loss mass.
            all_nl_mass_idx_start = self._generate_index(peak_data_nl["nl_mass"], max_indexed_mz)
            neutral_loss_index = (
                all_nl_mass_idx_start,
                peak_data_nl["nl_mass"],
                peak_data_nl["intensity"],
                peak_data_nl["spec_idx"],
                peak_data_nl["ion_mz"],  # The fragment ion m/z corresponding to the neutral loss mass.
            )
        else:
            neutral_loss_index = None

        self.index = fragment_ion_index, neutral_loss_index
        return self.index

    def get_topn_spec_idx_and_similarity(
        self,
        similarity_array,
        topn=None,
        min_similarity=0.1,
    ):
        """
        Retrieve the top-N spectra ranked by similarity score.

        Given an array of similarity scores, this method returns the indices and values of the best-matching spectra, optionally subject to minimum similarity filtering.

        Parameters
        ----------
        similarity_array : array-like
            A 1D array of similarity scores for all spectra in a group.

        topn : int or None, optional
            The maximum number of spectra to return.  
            If ``None`` (default), all spectra are considered.

        min_similarity : float or None, optional
            Minimum similarity threshold required for a match to be included.  
            Scores below this value are ignored.  
            If ``None``, all matches are accepted.  
            Default is ``0.1``.

        Returns
        -------
        tuple
            A tuple ``(result_idx, result)`` where:

            - ``result_idx`` is a list of spectrum indices sorted by decreasing
            similarity.
            - ``result`` is the corresponding list of similarity values.

            Only entries satisfying ``similarity >= min_similarity`` are returned.

        """

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
