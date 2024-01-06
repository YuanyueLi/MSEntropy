===========
Quick start
===========


Overview
========

.. code-block:: python

    from ms_entropy import FlashEntropySearch

    # Construct the FlashEntropySearch class
    entropy_search = FlashEntropySearch()

    # Step 1: Build the index from the library spectra
    entropy_search.build_index(spectral_library)

    # Step 2: Search the library
    entropy_similarity = entropy_search.search(
        precursor_mz = query_spectrum_precursor_mz, peaks = query_spectrum_peaks)

------------

In detail
=========

Suppose you have a spectral library, you need to format it like this:

.. code-block:: python

    import numpy as np
    spectral_library = [{
        "id": "Demo spectrum 1-A",
        "precursor_mz": 150.0,
        "peaks": [[100.0, 1.0], [101.0, 1.0], [103.0, 1.0]]
    }, {
        "id": "Demo spectrum 2-C",
        "precursor_mz": 250.0,
        "peaks": np.array([[200.0, 1.0], [101.0, 1.0], [202.0, 1.0]], dtype=np.float32),
        "XXX": "YYY",
    }, {
        "id": "Demo spectrum 3-B",
        "precursor_mz": 200.0,
        "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32),
        "metadata": "ABC"
    }, {
        "precursor_mz": 350.0,
        "peaks": [[100.0, 1.0], [101.0, 1.0], [302.0, 1.0]]}]

Note that the ``precursor_mz`` and ``peaks`` keys are required, the reset of the keys are optional.

Then you have your query spectrum looks like this:

.. code-block:: python

    query_spectrum = {"precursor_mz": 150.0,
                      "peaks": [[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]]}

You can call the ``FlashEntropySearch`` class to search the library like this:

.. code-block:: python

    from ms_entropy import FlashEntropySearch
    entropy_search = FlashEntropySearch()
    # Step 1: Build the index from the library spectra
    spectral_library_new = entropy_search.build_index(spectral_library)
    # Step 2: Search the library
    entropy_similarity = entropy_search.search(
        precursor_mz=150.0, peaks=[[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]])

.. warning::
    It is important to note that for efficient identity searching, all spectra in the ``spectral_library`` get **re-sorted** based on their precursor m/z values during the indexing process, and the ``search`` function returns the similarity scores in the same order as the **re-sorted spectra**. The ``build_index`` function returns a list of these **re-sorted spectra**, useful for mapping results back to the original spectra.

    In this example, the original order of the spectra is `"Demo spectrum 1-A", "Demo spectrum 2-C", "Demo spectrum 3-B", ...` After indexing, the spectra are re-sorted by precursor m/z, so the order becomes `"Demo spectrum 1-A", "Demo spectrum 3-B", "Demo spectrum 2-C", ...` You can check the difference by printing the ``spectral_library_new`` variable: ``print(spectral_library_new)``. The variable ``entropy_similarity`` will have the same order as ``spectral_library_new``.

After that, you can print the results like this:

.. code-block:: python

    import pprint
    pprint.pprint(entropy_similarity)

The result will look like this:

.. code-block:: python

    {'hybrid_search': array([0.6666666 , 0.99999994, 0.99999994, 0.99999994], dtype=float32),
    'identity_search': array([0.6666667, 0.       , 0.       , 0.       ], dtype=float32),
    'neutral_loss_search': array([0.6666666, 0.       , 0.6666666, 0.3333333], dtype=float32),
    'open_search': array([0.6666666 , 0.99999994, 0.3333333 , 0.6666666 ], dtype=float32)}

The values are the similarity scores for each spectrum in the ``spectral_library_new`` list. For example, the array ``[0.6666666 , 0.99999994, 0.3333333 , 0.6666666]`` in the ``open_search`` key means that the query spectrum has a similarity score of `0.6666666` with the first spectrum in the ``spectral_library_new`` list, which is **"Demo spectrum 1-A"** or ``entropy_search[0]``, a similarity score of `0.99999994` with the second spectrum in the ``spectral_library_new`` list, which is **"Demo spectrum 3-B"** or ``entropy_search[1]``, and a similarity score of `0.3333333` with the third spectrum in the ``spectral_library_new`` list (``entropy_search[2]``), and so on.

.. note::
    In default, the ``search`` function will return the similarity scores for all four search modes, which are ``identity_search``, ``open_search``, ``neutral_loss_search``, and ``hybrid_search``. To save time, you can specify the search mode by setting the ``method`` parameter, for example, ``method = {'identity', 'open'}`` will only return the similarity scores for ``identity_search`` and ``open_search``. `Click here <./entropy_search_basic_usage.html#search-function>`_ for more details.

------------

Examples
========

You can find several examples of how to use the package in the ``examples`` directory, the ``example.py`` script is a good starting point to get familiar with the package.

------------

Want more?
==========

Still have questions? Want more functions?

We also provided more function tools to help you calculate the spectral similarity, please go to the rest sections for more information.
