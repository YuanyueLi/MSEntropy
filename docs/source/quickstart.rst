===========
Quick start
===========


Super quick start
=================

Installation
------------
Python >= 3.7 is required.

.. code-block:: bash

    pip install ms_entropy


Basic usage for Flash entropy search
------------------------------------

.. code-block:: python

    from ms_entropy import FlashEntropySearch
    entropy_search = FlashEntropySearch()
    entropy_search.build_index(spectral_library)
    entropy_similarity = entropy_search.search(
        precursor_mz = query_spectrum_precursor_mz, peaks = query_spectrum_peaks)

------------

In detail
=========

Suppose you have a spectral library, you need to format it like this:

.. code-block:: python

    import numpy as np
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
    spectral_library = entropy_search.build_index(spectral_library)
    # Step 2: Search the library
    entropy_similarity = entropy_search.search(
        precursor_mz=150.0, peaks=[[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]])


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

------------

Examples
========

You can find several examples of how to use the package in the ``examples`` directory, the ``example.py`` script is a good starting point to get familiar with the package.

------------

Want more?
==========

Still have questions? Want more functions?

We also provided more function tools to help you calculate the spectral similarity, please go to the rest sections for more information.
