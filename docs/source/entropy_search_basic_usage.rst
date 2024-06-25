===========
Basic usage
===========

Flash Entropy Search is an efficient algorithm designed for quickly searching large spectral libraries. For details, see the reference: `Li, Y., Fiehn, O., Flash entropy search to query all mass spectral libraries in real time. 04 April 2023, PREPRINT (Version 1) available at Research Square. <https://doi.org/10.21203/rs.3.rs-2693233/v1>`_.


Step 0: Prepare the library spectra
===================================

Before you start, you need to prepare your library spectra. This data should be represented as a list of dictionaries, where each dictionary corresponds to a spectrum. Below is an example demonstrating the appropriate format for the library spectra:

.. code-block:: python

    import numpy as np
    spectral_library =  [{
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


In this format, each dictionary must contain the keys ``precursor_mz`` and ``peaks``. The ``precursor_mz`` key represents the precursor m/z, which should be a float. The ``peaks`` key represents the peaks of the spectrum, which can be a list of lists or a numpy array in the format ``[[mz, intensity], [mz, intensity], ...]``.

All other keys are optional and can be used to store any metadata associated with the spectrum. For instance, the ``id`` key could be used to identify the spectrum, while other keys could hold additional information.


Step 1: Building the index
==========================

The initial step in using the Flash Entropy Search algorithm involves building an index for the library spectra. The ``build_index`` function of the ``FlashEntropySearch`` class helps you accomplish this. This function takes your ``spectral_library`` as an argument. The spectra in the ``spectral_library`` do not require any pre-processing, as the index building process will handle this automatically.

.. code-block:: python

    from ms_entropy import FlashEntropySearch
    entropy_search = FlashEntropySearch()
    entropy_search.build_index(spectral_library)


.. warning::
    It is important to note that for efficient identity searching, all spectra in the ``spectral_library`` get **re-sorted** based on their precursor m/z values during the indexing process. The ``build_index`` function returns a list of these re-sorted spectra, useful for mapping results back to their original order.

    If you need to maintain the original order of spectra, consider adding metadata, such as an ``id`` field, to your spectra. This metadata will help keep track of the original sequence even after re-sorting.


.. note::
    If you want to calcualte unweighted entropy similarity instead of weighted entropy similarity, you can set the ``intensity_weight`` parameter to ``None`` when constructing the ``FlashEntropySearch`` object. For example:

    .. code-block:: python

        entropy_search = FlashEntropySearch(intensity_weight=None)
        entropy_search.build_index(spectral_library)
        ...


Step 2: Searching the library
=============================

After building the index, you can search your query spectrum against the library spectra using the ``search`` function. This function takes the query spectrum as input, returning the similarity score of each spectrum in the library, in the same order as the spectral library returned by the ``build_index`` function.

.. _search-function:

The ``search`` function accepts the following parameters:

- ``precursor_mz``: The precursor m/z of the query spectrum.

- ``peaks``: The peaks of the query spectrum, which can be either a list of lists or a numpy array, which in format of ``[[mz, intensity], [mz, intensity], ...]``.

- ``ms1_tolerance_in_da``: The mass tolerance to apply to the precursor m/z in Da, used only for identity search.

- ``ms2_tolerance_in_da``: The mass tolerance to apply to the fragment ions in Da.

- ``method``: The search method to employ. Available methods include ``identity``, ``open``, ``neutral_loss``, and ``hybrid``. You can use a string or a list/set of these four strings like ``{'identity', 'open'}``. Use ``all`` to apply all four methods. The default value is ``all``.

- ``target``: Where to perform the similarity calculation, on the ``cpu`` or ``gpu``. The default value is ``cpu``.

- ``precursor_ions_removal_da``: The mass tolerance for removing the precursor ions in Da. Fragment ions with m/z larger than ``precursor_mz - precursor_ions_removal_da`` will be removed. Based on our tests, removing precursor ions can enhance search performance.

- ``noise_threshold``: The intensity threshold for removing noise peaks. Peaks with intensity smaller than ``noise_threshold * max(fragment ion's intensity)`` will be removed.

- ``max_num_peaks``: Keep only the top ``max_num_peaks`` peaks with the highest intensity.

To run the search function, use the following code:

.. code-block:: python

    entropy_similarity = entropy_search.search(
        precursor_mz = 150.0,
        peaks = [[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]]
    )


The ``search`` function returns a dictionary, where the key is the search method and the value is a list of similarity scores. The scores align with the order of the spectral library returned by the ``build_index`` function. An example of the results is shown below:

.. code-block:: python

    {
        'identity_search': [0.0, 0.5, 0.0, 0.8],
        'open_search': [0.0, 0.0, 0.3, 0.8],
        'neutral_loss_search': [0.2, 0.0, 0.7, 0.0],
        'hybrid_search': [0.2, 0.5, 1.0, 0.8]
    }


Alternative: individual search functions
========================================

Instead of using the ``search`` function that automatically (1) cleans the query spectrum and (2) performs the library search, you have the option to manually perform these tasks in two separate steps. You can use the ``clean_spectrum_for_search`` function to clean the query spectrum and then use individual search functions to search the library. Both approaches are equivalent, and you can choose the one that suits you best.

Clean the query spectrum
------------------------

Before performing a library search, the query spectrum should be pre-processed using the ``clean_spectrum_for_search`` function. This function accomplishes the following:

1. Remove empty peaks (m/z <= 0 or intensity <= 0).

2. Remove peaks with m/z values greater than ``precursor_mz - precursor_ions_removal_da`` (removes precursor ions to improve the quality of spectral comparison).

3. Centroid the spectrum by merging peaks within +/- ``min_ms2_difference_in_da`` and sort the resulting spectrum by m/z.

4. Remove peaks with intensity less than ``noise_threshold`` * maximum intensity.

5. Retain only the top max_peak_num peaks and remove all others.

6. Normalize the intensity to sum to 1.

Assuming you have your query spectrum as:

.. code-block:: python

    query_spectrum = {"precursor_mz": 150.0,
                      "peaks": [[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]]}

To utilize the ``clean_spectrum_for_search`` function, call it on your query spectrum, passing in the relevant parameters:

.. code-block:: python

    query_spectrum['peaks'] = entropy_search.clean_spectrum_for_search(
        precursor_mz = query_spectrum['precursor_mz'],
        peaks = query_spectrum['peaks']
    )

We also provide a separate function called ``clean_spectrum`` that performs the same cleaning steps as ``clean_spectrum_for_search``. Here's how to call this function:

.. code-block:: python

    from ms_entropy import clean_spectrum
    precursor_ions_removal_da = 1.6
    query_spectrum['peaks'] = clean_spectrum(
        peaks = query_spectrum['peaks'],
        max_mz = query_spectrum['precursor_mz'] - precursor_ions_removal_da
    )

Both of these functions serve the same purpose and can be used interchangeably. You can select the one that suits your needs.


Performing library search using individual search functions
-----------------------------------------------------------

There are four individual search functions available for library searching:

- ``identity_search`` for Identity search
- ``open_search`` for Open search
- ``neutral_loss_search`` for Neutral loss search
- ``hybrid_search`` for Hybrid search

Each of these functions takes the ``pre-cleaned`` query spectrum as input, along with the spectral library index built in Step 1, and returns the similarity score for each spectrum in the library, in the same order as the spectral library that was returned by the ``set_library_spectra`` function.

.. warning::
    Remember, when using any of the four individual search functions, the peaks must be pre-processed by either the ``clean_spectrum_for_search`` or the ``clean_spectrum`` function. Failure to do so will result in an error.

Each search function accepts the following parameters:

- ``precursor_mz``: The precursor m/z of the query spectrum.
- ``peaks``: The peaks of the query spectrum.
- ``ms1_tolerance_in_da``: The mass tolerance to use for the precursor m/z in Da.
- ``ms2_tolerance_in_da``: The mass tolerance to use for the fragment ions in Da.
- ``target``: Specifies whether to run the similarity calculation on CPU or GPU. The default value is ``cpu``.

Here's an example of how you can use these functions:

.. code-block:: python
        
    # Identity search
    entropy_similarity = entropy_search.identity_search(
        precursor_mz = query_spectrum['precursor_mz'],
        peaks = query_spectrum['peaks'],
        ms1_tolerance_in_da = 0.01,
        ms2_tolerance_in_da = 0.02
    )

    # Open search
    entropy_similarity = entropy_search.open_search(
        peaks = query_spectrum['peaks'],
        ms2_tolerance_in_da = 0.02
    )

    # Neutral loss search
    entropy_similarity = entropy_search.neutral_loss_search(
        precursor_mz = query_spectrum['precursor_mz'],
        peaks = query_spectrum['peaks'],
        ms2_tolerance_in_da = 0.02
    )

    # Hybrid search
    entropy_similarity = entropy_search.hybrid_search(
        precursor_mz = query_spectrum['precursor_mz'],
        peaks = query_spectrum['peaks'],
        ms2_tolerance_in_da = 0.02
    )
