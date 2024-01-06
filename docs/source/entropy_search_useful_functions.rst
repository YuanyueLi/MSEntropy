================
Useful functions
================

Reading Spectra from a File
===========================
For ease of use, a function named ``read_one_spectrum`` is provided in the ``ms_entropy`` package, allowing you to easily read spectra from a file. Here is an example of how you can use it:

.. code-block:: python

    from ms_entropy import read_one_spectrum
    for spectrum in read_one_spectrum('path/to/spectrum/file'):
        print(spectrum)

This function returns a dictionary, where each key-value pair corresponds to a specific metadata of the spectrum.

Currently, the ``read_one_spectrum`` function supports the following file formats: ``.mgf``, ``.msp``, ``.mzML``, and ``.lbm2`` from the MS-DIAL software.

----------------

Get the top-n results from the Flash entropy search results
===========================================================

Once you have conducted a search in your spectral library, you may want to focus only on the top-N results, or the results with a similarity score that is higher than a certain threshold. The ``get_topn_matches`` function has been designed specifically for this purpose.

The ``get_topn_matches`` function takes three parameters:

- ``similarity_array``: The similarity scores that the search function has returned.

- ``topn``: The number of top results you want to retrieve. If you set this to None, all results will be retrieved.

- ``min_similarity``: The minimum similarity score that results should have. If you set this to None, all results will be retrieved.

The function will return a list of dictionaries. Each dictionary corresponds to a spectrum in the library. The dictionary is similar to the one in the library spectra (the input of the ``build_index``), with the addition of an ``entropy_similarity`` key to store the similarity score of the spectrum.

Here's how you can use the ``get_topn_matches`` function:

.. code-block:: python

    topn_match = entropy_search.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)


This example will return a list of the top 3 matches with a similarity score greater than 0.01.

----------------

Get the metadata of a specifical spectrum from the Flash entropy search object
==============================================================================

After you've conducted a search in your spectral library, you may want to retrieve the metadata of a specific spectrum. For this, you can use the ``__getitem__`` function.

For instance, let's say that after a search, you found that the third spectrum (index start from 0) in the library has the highest similarity score. You can call ``entropy_search[2]`` to retrieve the metadata of the third spectrum.

Here's an example of how you can use the ``__getitem__`` function:

.. code-block:: python

    from ms_entropy import FlashEntropySearch
    entropy_search = FlashEntropySearch()
    entropy_search.build_index(spectral_library)

    # Get the metadata of the third spectrum
    metadata = entropy_search[2]

The metadata was extracted and stored when you called the ``build_index`` function. The data will remain available even if you save and reload the index using either the pickle module or the read and write functions.

----------------

Get the matched peaks number of query spectrum to the library Spectra
=====================================================================

If you also want to know the number of matched peaks between the query spectrum and the library spectra, you can set the ``get_matched_peaks_number`` parameters to ``True``. Then, the returned results will be a list of two numpy arrays. The first array contains the similarity scores, and the second array contains the number of matched peaks.

At this moment, the ``get_matched_peaks_number`` parameter is only supported by the ``identity_search``, ``open_search``, and ``neutral_loss_search`` functions.

Here's an example of how you can use the ``get_matched_peaks_number`` parameter:

.. code-block:: python

    import numpy as np
    from ms_entropy import FlashEntropySearch

    spectral_library = [{
        "id": "Demo spectrum 1",
        "precursor_mz": 150.0,
        "peaks": [[100.0, 1.0], [101.0, 1.0], [103.0, 1.0]]
    }, {
        "id": "Demo spectrum 2",
        "precursor_mz": 200.0,
        "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32),
        "metadata": "ABC"
    }, {
        "id": "Demo spectrum 3",
        "precursor_mz": 250.0,
        "peaks": np.array([[200.0, 1.0], [101.0, 1.0], [202.0, 1.0]], dtype=np.float32),
        "XXX": "YYY",
    }, {
        "precursor_mz": 350.0,
        "peaks": [[100.0, 1.0], [101.0, 1.0], [302.0, 1.0]]}]
    query_spectrum = {"precursor_mz": 150.0,
                      "peaks": [[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]]}

    entropy_search = FlashEntropySearch()
    # Step 1: Build the index from the library spectra
    spectral_library = entropy_search.build_index(spectral_library)
    # Step 2: Clean the query spectrum
    query_spectrum['peaks'] = entropy_search.clean_spectrum_for_search(
        precursor_mz = query_spectrum['precursor_mz'],
        peaks = query_spectrum['peaks']
    )
    # Step 3: Search the library
    # This parameter is supported by the identity_search, open_search, and neutral_loss_search functions
    entropy_similarity, matched_peaks_number = entropy_search.identity_search(
        precursor_mz = query_spectrum['precursor_mz'],
        peaks = query_spectrum['peaks'],
        ms1_tolerance_in_da = 0.01,
        ms2_tolerance_in_da = 0.02,
        output_matched_peak_number = True
    )
    print(entropy_similarity)
    print(matched_peaks_number)

----------------

Save and load index for the Flash entropy search object
=======================================================

After you have built the index, you have the option to save it to disk for later use.

Using pickle
------------

You can use Python's built-in ``pickle`` module to save and load the ``FlashEntropySearch`` object, as follows:

.. code-block:: python

    import pickle
    # Save the index
    with open('path/to/index', 'wb') as f:
        pickle.dump(entropy_search, f)
    # And load the index
    with open('path/to/index', 'rb') as f:
        entropy_search = pickle.load(f)

Using ``read`` and ``write`` functions
--------------------------------------

We also provide ``read`` and ``write`` functions to save and load the index.


To save a ``FlashEntropySearch`` object to disk:

.. code-block:: python

    entropy_search.write('path/to/index')


To load a ``FlashEntropySearch`` object from disk:

.. code-block:: python

    entropy_search = FlashEntropySearch()
    entropy_search.read('path/to/index')


If you're working with a very large spectral library, or your computer's memory is limited, you can use the ``low_memory`` parameter to partially load the library and reduce the memory usage. For example:

.. code-block:: python

    entropy_search = FlashEntropySearch(low_memory=True)
    entropy_search.read('path/to/index')

The index only needs to be built once. After that, you can use the read function to load the index. If you built the index using the ``low_memory=False`` mode, you can still load it using a ``FlashEntropySearch`` object with either the ``low_memory=False`` or ``low_memory=True`` mode.
