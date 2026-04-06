===========
Quick start
===========
Examples of DynamicEntropySearch can be found in `Dynamic Entropy Search Github Repository <https://github.com/2bereal-me/DynamicEntropySearch>`_.

.. code-block:: python

    from ms_entropy import DynamicEntropySearch

    # Step 1: Construct the DynamicEntropySearch class.
    entropy_search = DynamicEntropySearch(path_data=path_of_your_library)

    # Step 2: Construct or update the index from the library spectra.
    entropy_search.add_new_spectra(spectra_list=spectra_1_for_library)
    entropy_search.add_new_spectra(spectra_list=spectra_2_for_library)
    ......

    # Step 3: Call build_index() and write() lastly to end adding operation.
    entropy_search.build_index()
    entropy_search.write()

    # Step 4: Perform search.
    entropy_similarity=entropy_search.search(
        precursor_mz=query_spectrum_precursor_mz,
        peaks=query_spectrum_peaks)



