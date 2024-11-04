========
Examples
========

This package provides several examples in the Github root directory to help users understand its usage. These examples cover a variety of use cases, including building an index, searching spectra, and evaluating the performance of search functions.

These examples can be found in the ``examples`` directory. Here is a brief description of some of the available examples:


example.py
==========

This script demonstrates how to use the Flash Entropy Search module from scratch. It showcases the process of performing similarity searches using various methods. The output should look like the following:

.. code-block:: bash

    -------------------- All types of similarity search --------------------
    {'hybrid_search': array([0.6666666 , 0.99999994, 0.99999994, 0.99999994], dtype=float32),
    'identity_search': array([0.6666667, 0.       , 0.       , 0.       ], dtype=float32),
    'neutral_loss_search': array([0.6666666, 0.       , 0.6666666, 0.3333333], dtype=float32),
    'open_search': array([0.6666666 , 0.99999994, 0.3333333 , 0.6666666 ], dtype=float32)}


example-2.py
============

The ``example-2.py`` script provides a detailed illustration of how to utilize the Flash entropy search module from scratch. Special focus is placed on dealing with metadata and retrieving the top-n results. The script's output will resemble the following:

.. code-block:: bash

    -------------------- Identity search --------------------
    [{'entropy_similarity': 0.6666667,
    'id': 'Demo spectrum 1',
    'precursor_mz': 150.0}]
    -------------------- Open search --------------------
    [{'entropy_similarity': 0.99999994,
    'id': 'Demo spectrum 2',
    'metadata': 'ABC',
    'precursor_mz': 200.0},
    {'entropy_similarity': 0.6666666, 'precursor_mz': 350.0},
    {'entropy_similarity': 0.6666666,
    'id': 'Demo spectrum 1',
    'precursor_mz': 150.0}]
    -------------------- Neutral loss search --------------------
    [{'XXX': 'YYY',
    'entropy_similarity': 0.6666666,
    'id': 'Demo spectrum 3',
    'precursor_mz': 250.0},
    {'entropy_similarity': 0.6666666,
    'id': 'Demo spectrum 1',
    'precursor_mz': 150.0},
    {'entropy_similarity': 0.3333333, 'precursor_mz': 350.0}]
    -------------------- Hybrid search --------------------
    [{'entropy_similarity': 0.99999994, 'precursor_mz': 350.0},
    {'XXX': 'YYY',
    'entropy_similarity': 0.99999994,
    'id': 'Demo spectrum 3',
    'precursor_mz': 250.0},
    {'entropy_similarity': 0.99999994,
    'id': 'Demo spectrum 2',
    'metadata': 'ABC',
    'precursor_mz': 200.0}]


example-with_individual_search_function.py
==========================================

This example demonstrates how to use the Flash entropy search from scratch, but unlike ``example-2.py``, it leverages individual search functions rather than the composite ``search`` function. The expected output should be identical to ``example-2.py``.


example-search_mona-with_read_write_functions.py
================================================

This example demonstrates the usage of the Flash entropy search to perform comprehensive searches of the entire [MassBank.us (MoNA)](https://massbank.us/) database.

Upon its initial execution, the script may take approximately 10-20 minutes to download and parse spectra from the MoNA .msp files. An additional 2-4 minutes will be required to build the index for the MoNA library. However, during subsequent runs, the index will be fetched directly from disk, significantly reducing the runtime to less than a second. Here's what you can expect the output to look like:

.. code-block:: bash

    Loading index
    Downloading https://mona.fiehnlab.ucdavis.edu/rest/downloads/retrieve/03d5a22c-c1e1-4101-ac70-9a4eae437ef5 to /p/FastEntropySearch/github_test/data/mona-2023-03-23.zip
    Loading spectra from /p/FastEntropySearch/github_test/data/mona-2023-03-23.zip, this may take a while.
    Loaded 811840 positive spectra and 1198329 negative spectra.
    Building index, this will only need to be done once.
    Building index for spectra with ion mode P
    Building index for spectra with ion mode N
    Saving index
    ********************************************************************************
    Identity Search with Flash Entropy Search
    Finished identity search in 0.0017 seconds with 1196680 results.
    Top 5 matches:
    Rank 1: AU116754 with score 1.0000
    Rank 2: AU116755 with score 0.8081
    Rank 3: AU116753 with score 0.6565
    Rank 4: AU116752 with score 0.2717
    ********************************************************************************
    Open Search with Flash Entropy Search
    Finished open search in 0.0006 seconds with 1196680 results.
    Top 5 matches:
    Rank 1: AU116754 with score 1.0000
    Rank 2: AU116755 with score 0.8081
    Rank 3: AU116753 with score 0.6565
    Rank 4: CCMSLIB00004751228 with score 0.4741
    Rank 5: LU040151 with score 0.4317
    ********************************************************************************
    Neutral Loss Search with Flash Entropy Search
    Finished neutral loss search in 0.0006 seconds with 1196680 results.
    Top 5 matches:
    Rank 1: AU116754 with score 1.0000
    Rank 2: AU116755 with score 0.8081
    Rank 3: AU116753 with score 0.6565
    Rank 4: LipidBlast2022_1230911 with score 0.3796
    Rank 5: LipidBlast2022_1230977 with score 0.3796
    ********************************************************************************
    Hybrid Search with Flash Entropy Search
    Finished hybrid search in 0.0010 seconds with 1196680 results.
    Top 5 matches:
    Rank 1: AU116754 with score 1.0000
    Rank 2: AU116755 with score 0.8081
    Rank 3: AU116753 with score 0.6565
    Rank 4: CCMSLIB00004751228 with score 0.4741
    Rank 5: LU040151 with score 0.4317


example-search_mona-with_read_write_functions-low_memory_usage.py
=================================================================

This example, demonstrates the usage of the Flash entropy search to search the [MassBank.us (MoNA)](https://massbank.us/) database, utilizing less memory than the ``example-search_mona-with_read_write_functions.py`` script.

To get started, either ``example-search_mona-with_read_write_functions.py`` or ``example-search_mona-with_read_write_functions-low_memory_usage.py`` should be executed to build the index for the MoNA library. On subsequent runs of either scripts, the previously constructed index will be fetched directly from disk, resulting in faster load times.

Upon successful index construction, a second run of ``example-search_mona-with_read_write_functions.py`` on my computer consumed about 1,212MB memory to search a single spectrum against the entire MassBank.us (MoNA) library. Conversely, running ``example-search_mona-with_read_write_functions-low_memory_usage.py`` for the second time only required about 84MB of memory to perform the same task. This lower memory usage is especially beneficial when dealing with an extensive spectral library and you have limited computer memory.


example-search_mona-with_pickle_functions.py
============================================

An other example shows how to use the Flash entropy search to search the [MassBank.us (MoNA)](https://massbank.us/) database. This example use built-in pickle functions to save and load index.
