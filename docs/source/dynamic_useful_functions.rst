================
Useful functions
================

Here we introduce functions that you may find useful when using Dynamic Entropy Search.


About index construction
==========================

Manually convert the index to a Flash Entropy Search-compatible format
-----------------------------------------------------------------------
Based on the Flash Entropy Search, the index of Dynamic Entropy Search can be converted to a compact structure which can be used in Flash Entropy Search, resulting in a faster search speed and less storage.

.. code-block:: python

    from ms_entropy import DynamicEntropySearch
    from pathlib import Path
    # Choose an existing index and assign the path
    entropy_search = DynamicEntropySearch(path_data=path_of_your_library)
    
    # Manually sort the blocks in index
    entropy_search.convert_to_fast_search()

    # Manually convert the index to a compact structure
    entropy_search.read()
    group_num=len(entropy_search.group_start)
    for i in range(group_num):
        group_path=Path(path_of_your_library)/f"{i}"
        entropy_search.entropy_search.path_data=group_path
        entropy_search.entropy_search.read()
        entropy_search.convert_current_index_to_flash()

This operation internally sorts the blocks in the index and then removes reserved space in the index. 

After this process, index can be converted to the structure compatible with Flash Entropy Search. Performance of search will improve when using the search functions in ``DynamicEntropySearch``.

.. note::
    When using ``add_new_spectra()`` functions to update index, there will be an automatic conversion of index structure if the size of this group meets the limit.
    You can set ``convert_to_flash`` in ``add_new_spectra()`` and ``build_index()`` to ``False`` to disable this feature.
    
    .. code-block:: python

        from ms_entropy import DynamicEntropySearch
        entropy_search = DynamicEntropySearch(path_data=path_of_your_library)
        entropy_search.add_new_spectra(spectra_list=spectra_1_for_library, convert_to_flash=False)
        entropy_search.add_new_spectra(spectra_list=spectra_2_for_library, convert_to_flash=False)
        entropy_search.add_new_spectra(spectra_list=spectra_3_for_library, convert_to_flash=False)
        ......
        entropy_search.build_index(convert_to_flash=False)
        entropy_search.write()

    It should be noted that if the index of group is already a compact structure, converting operation and subsequent adding operation will both result in an error.

Construct an index only for open search
-----------------------------------------

If you only need to construct index for open search, it is unnecessary to process neutral loss data.
By setting ``index_for_neutral_loss`` in ``add_new_spectra()`` and ``build_index()`` to ``False``, you can construct the index for open search more efficiently.

Here's an example:

.. code-block:: python

    from ms_entropy import DynamicEntropySearch
    entropy_search = DynamicEntropySearch(path_data=path_of_your_library)
    entropy_search.add_new_spectra(spectra_list=spectra_1_for_library, index_for_neutral_loss=False)
    entropy_search.add_new_spectra(spectra_list=spectra_2_for_library, index_for_neutral_loss=False)
    entropy_search.add_new_spectra(spectra_list=spectra_3_for_library, index_for_neutral_loss=False)
    ......
    entropy_search.build_index(index_for_neutral_loss=False)
    entropy_search.write()

.. warning::
    Once ``index_for_neutral_loss`` is set to ``False``, it will no longer be possible to construct neutral loss index of this library. Keep this parameter to ``False`` all the time to avoid errors.
    What's more, only open search can be performed under this circumstance. It is necessary to check the value of ``method`` when using search functions. Performing identity search, neutral loss search or hybrid search can cause error.


About spectra search
==========================


Different search functions can serve different objectives.

Make sure the query spectra are all cleaned before search. 


General search
---------------

Search top matches with metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from ms_entropy import DynamicEntropySearch

    # Assign the path for your library
    entropy_search=DynamicEntropySearch(path_data=path_of_your_library)

    # For query_spectra_list, iterate it to perform search for each elements.
    for spec in query_spectra_list:
        result=entropy_search.search_topn_matches(
                precursor_mz=spec['precursor_mz'],
                peaks=spec['peaks'],
                ms1_tolerance_in_da=0.01, # You can change ms1_tolerance_in_da as needed.
                ms2_tolerance_in_da=0.02, # You can change ms2_tolerance_in_da as needed.
                method='open', # or 'neutral_loss' or 'hybrid' or 'identity'.
                clean=True, # If you don't want to use the internal clean process in this function, set it to False.
                topn=3, # You can change topn as needed.
                need_metadata=True, # Set it to True if need metadata.
        )
        # After that, you can print the result like this:
        print(result)

An example result:

.. code-block:: bash

    [{
    'id': 'Demo spectrum 2', 
    'precursor_mz': 200.0, 
    'peaks': array([[100.        ,   0.33333334], [101.        ,   0.33333334], [102.        ,   0.33333334]], dtype=float32), 
    'metadata': 'ABC', 
    'open_search_entropy_similarity': np.float32(0.99999994)
    }, {
    'id': 'Demo spectrum 1', 
    'precursor_mz': 150.0, 
    'peaks': array([[100.        ,   0.33333334], [101.        ,   0.33333334], [103.        ,   0.33333334]], dtype=float32), 
    'open_search_entropy_similarity': np.float32(0.6666666)
    }, {
    'precursor_mz': 350.0, 
    'peaks': array([[100.        ,   0.33333334], [101.        ,   0.33333334],[302.        ,   0.33333334]], dtype=float32), 'open_search_entropy_similarity': np.float32(0.6666666)
    }]



Search that only requires similarity without metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from ms_entropy import DynamicEntropySearch

    # Assign the path for your library
    entropy_search=DynamicEntropySearch(path_data=path_of_your_library)

    ### Use `search()` and get an array with all entropy similarities ###
    result=entropy_search.search(
            precursor_mz=query_spectrum['precursor_mz'],
            peaks=query_spectrum['peaks'],
            ms1_tolerance_in_da=0.01, # You can change ms1_tolerance_in_da as needed.
            ms2_tolerance_in_da=0.02, # You can change ms2_tolerance_in_da as needed.
            method='all', # or 'neutral_loss' or 'hybrid' or 'identity' or 'open'.
            clean=True, # If you don't want to use the internal clean process in this function, set it to False.
    )
    print(result)

Example result:

.. code-block:: bash

    {
    'identity_search': array([0.6666666, 0.       , 0.       , 0.       , 0.       , 0.       ], dtype=float32), 
    'open_search': array([0.6666666 , 0.99999994, 0.3333333 , 0.6666666 , 0.        , 0.        ], dtype=float32), 
    'neutral_loss_search': array([0.6666666, 0.       , 0.6666666, 0.3333333, 0.       , 0.       ], dtype=float32), 
    'hybrid_search': array([0.6666666 , 0.99999994, 0.99999994, 0.99999994, 0.        , 0.        ], dtype=float32)
    }

This result:

- includes the results of all search methods because method is set to all.
- returns only similarity array in the order of spectra in the library.



Specific search
--------------------
Specific searches include identity search, open search, neutral loss search and hybrid search.

.. note::
    These following specific search functions don't have an internal clean function. 
    Thus, external clean function is necessary before using these search functions. See :ref:`external clean before search <clean-spectrum-function>` and the following examples.

Identity search
^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from ms_entropy import DynamicEntropySearch

    # Assign the path for your library
    entropy_search=DynamicEntropySearch(path_data=path_of_your_library)

    precursor_ions_removal_da = 1.6

    query_spectrum['peaks'] = clean_spectrum(
        peaks = query_spectrum['peaks'],
        max_mz = query_spectrum['precursor_mz'] - precursor_ions_removal_da
    )

    ### Use `identity_search()` and get an array with all entropy similarities based on identity search ###
    result=entropy_search.identity_search(
            precursor_mz=query_spectrum['precursor_mz'],
            peaks=query_spectrum['peaks'],
            ms1_tolerance_in_da=0.01, # You can change ms1_tolerance_in_da as needed.
            ms2_tolerance_in_da=0.02, # You can change ms2_tolerance_in_da as needed.
    )
    print(result)

Example result:

.. code-block:: bash

    [0.         0.         0.99999994 0.         0.         0.        ]

This result:

- includes the results of identity search.
- returns only similarity array in the order of spectra in the library.


Open search
^^^^^^^^^^^^^^^^^^
.. code-block:: python

    from ms_entropy import DynamicEntropySearch

    # Assign the path for your library
    entropy_search=DynamicEntropySearch(path_data=path_of_your_library)

    query_spectrum['peaks'] = clean_spectrum(
        peaks = query_spectrum['peaks'],
    )

    ### Use `open_search()` and get an array with all entropy similarities based on open search ###
    result=entropy_search.open_search(
            peaks=query_spectrum['peaks'],
            ms2_tolerance_in_da=0.02, # You can change ms2_tolerance_in_da as needed.
    )
    print(result)


Example result:

.. code-block:: bash

    [0.3333333  0.3333333  0.99999994 0.3333333  0.         0.        ]


This result:

- includes the results of open search.
- returns only similarity array in the order of spectra in the library.


Neutral loss search
^^^^^^^^^^^^^^^^^^
.. code-block:: python

    from ms_entropy import DynamicEntropySearch

    # Assign the path for your library
    entropy_search=DynamicEntropySearch(path_data=path_of_your_library)

    precursor_ions_removal_da = 1.6

    query_spectrum['peaks'] = clean_spectrum(
        peaks = query_spectrum['peaks'],
        max_mz = query_spectrum['precursor_mz'] - precursor_ions_removal_da
    )

    ### Use `neutral_loss_search()` and get an array with all entropy similarities based on neutral loss search ###
    result=entropy_search.neutral_loss_search(
            precursor_mz=query_spectrum['precursor_mz'],
            peaks=query_spectrum['peaks'],
            ms2_tolerance_in_da=0.02, # You can change ms2_tolerance_in_da as needed.
    )
    print(result)


Example result:

.. code-block:: bash

    [0.3333333  0.         0.99999994 0.3333333  0.         0.        ]


This result:

- includes the results of neutral loss search.
- returns only similarity array in the order of spectra in the library.

Hybrid search
^^^^^^^^^^^^^^^^^^
.. code-block:: python

    from ms_entropy import DynamicEntropySearch

    # Assign the path for your library
    entropy_search=DynamicEntropySearch(path_data=path_of_your_library)

    precursor_ions_removal_da = 1.6

    query_spectrum['peaks'] = clean_spectrum(
        peaks = query_spectrum['peaks'],
        max_mz = query_spectrum['precursor_mz'] - precursor_ions_removal_da
    )

    ### Use `hybrid_search()` and get an array with all entropy similarities based on hybrid search ###
    result=entropy_search.hybrid_search(
            precursor_mz=query_spectrum['precursor_mz'],
            peaks=query_spectrum['peaks'],
            ms2_tolerance_in_da=0.02, # You can change ms2_tolerance_in_da as needed.
    )
    print(result)

Example result:

.. code-block:: bash

    [0.6666666  0.3333333  0.99999994 0.6666666  0.         0.        ]

This result:

- includes the results of hybrid search.
- returns only similarity array in the order of spectra in the library.
