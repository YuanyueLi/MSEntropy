=================================
Basic usage - Search
=================================

Spectra comparison is performed as searching a query spectrum against the index in the reference library. 
You can perform identity search, open search, neutral loss search or hybrid search based on your need. 

Search with internal clean function
===================================

Suppose you have established a library locally under ``path_of_your_library`` using the aforementioned method.

Now you can perform search with a query spectrum in correct format like this:

.. code-block:: python

    import numpy as np
    # For each query spectrum, 'precursor_mz' and 'peaks' are necessary. 
    # 'precursor_mz' should be a float, and 'peaks' should be a 2D np.ndarray like np.ndarray([[m/z, intensity], [m/z, intensity], [m/z, intensity]...], dtype=np.float32).

    query_spectrum = {"precursor_mz": 150.0,
                    "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32)}


If your query spectra is a list consisting of several spectra:

.. code-block:: python

    import numpy as np

    # For each query_spectra_list, it is a list consisting of multiple dictionaries of query MS2 spectra.

    # For each query spectrum, 'precursor_mz' and 'peaks' are necessary. 
    # 'precursor_mz' should be a float, and 'peaks' should be a 2D np.ndarray like np.ndarray([[m/z, intensity], [m/z, intensity], [m/z, intensity]...], dtype=np.float32).

    query_spectra_list = [{
                    "precursor_mz": 150.0,
                    "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32)
                    },{
                    "precursor_mz": 250.0,
                    "peaks": np.array([[108.0, 1.0], [113.0, 1.0], [157.0, 1.0]], dtype=np.float32)
                    },{
                    "precursor_mz": 299.0,
                    "peaks": np.array([[119.0, 1.0], [145.0, 1.0], [157.0, 1.0]], dtype=np.float32)
                    },
                    ]

You can call the ``DynamicEntropySearch`` class with corresponding ``path_data`` to search the library like this:

.. code-block:: python

    from ms_entropy import DynamicEntropySearch

    # Select the path for your library
    entropy_search=DynamicEntropySearch(path_data=path_of_your_library)

    # Search the library and you can fetch the metadata from the results with the highest scores
    result=entropy_search.search_topn_matches(
            precursor_mz=query_spectrum['precursor_mz'],
            peaks=query_spectrum['peaks'],
            ms1_tolerance_in_da=0.01, # You can change ms1_tolerance_in_da as needed.
            ms2_tolerance_in_da=0.02, # You can change ms2_tolerance_in_da as needed.
            method='open', # Or 'neutral_loss' or 'hybrid' or 'identity'.
            precursor_ions_removal_da=1.6, # Peaks with m/z greater than ``precursor_mz - precursor_ions_removal_da`` are removed during cleaning. 
            noise_threshold=0.01, # Relative intensity threshold for noise filtering during cleaning. Peaks with intensity ``< noise_threshold * max(intensity)`` are removed.  
            min_ms2_difference_in_da=0.05, # Minimum spacing allowed between MS/MS peaks during cleaning.  
            max_peak_num=None, # Maximum number of peaks to keep after cleaning.  
            clean=True, # If you don't want to use the internal clean process in this function, set it to False.
            topn=3, # You can change topn as needed.
            need_metadata=True, # Set it to True if need metadata.
    )

    # After that, you can print the result like this:
    print(result)

.. note::

    - Cleaning the query spectrum is necessary. You can use the internal clean function of ``search_topn_matches()`` or seperate the clean and search process. This is introduced in the following part.

    - ``search_topn_matches()`` is suitable for identification that requires metadata when ``need_metadata`` in it is ``True``. If it is set to ``False``, the location of matched spectra in library will be returned.


The ``search_topn_matches()`` function accepts the following parameters:

- ``peaks``: The peaks of the query spectrum, which is a numpy array in format of ``[[mz, intensity], [mz, intensity], ...]``.

- ``precursor_mz``: The precursor m/z of the query spectrum.

- ``ms1_tolerance_in_da``: The mass tolerance to apply to the precursor m/z in Da, used only for identity search. Default is ``0.01``.

- ``ms2_tolerance_in_da``: The mass tolerance to apply to the fragment ions in Da. Default is ``0.02``.

- ``method``: The search method to employ. Available methods include ``identity``, ``open``, ``neutral_loss``, and ``hybrid``. A string is acceptable. The default value is ``open``.

- ``clean``: Whether to clean the query spectrum before searching. Default is ``True``.

- ``precursor_ions_removal_da``: The mass tolerance for removing the precursor ions in Da. Fragment ions with m/z larger than ``precursor_mz - precursor_ions_removal_da`` will be removed. Based on our tests, removing precursor ions can enhance search performance. Default is ``1.6``.

- ``noise_threshold``: The intensity threshold for removing noise peaks. Peaks with intensity smaller than ``noise_threshold * max(fragment ion's intensity)`` will be removed. Default is ``0.01``.

- ``min_ms2_difference_in_da``:  Minimum spacing allowed between MS/MS peaks during cleaning. Default is ``0.05`` Da.

- ``max_peak_num``: Maximum number of peaks to keep after cleaning. ``None`` keeps all peaks. Default is ``None``.

- ``topn``: Number of top-matching spectra to return. If ``None``, all spectra are returned. Default is ``3``.

- ``need_metadata``: If ``True`` (default), return the metadata dictionary for each matched spectrum. If ``False``, return `(global_index, similarity)` tuples instead. Default is ``True``.


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

In this result:

- This is generated from searching the query spectrum against an existing library. Select correct ``method`` to perform search.

- 3 top matched spectra are given in descending order of similarity. This is set by ``topn`` in ``search_topn_matches()``. If number of spectra with similarity greater than 0 is less than ``topn``, then output the actual matching number of spectra.

- Metadata of spectra is given if ``need_metadata`` in ``search_topn_matches()`` is set to ``True``. Users can add information for spectra when constructing library. These additional information other than 'precursor_mz' and 'peaks', like 'id', benefit the compound identification.

If the query spectra is a list, iterate it to perform search.

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


.. _external-clean-before-search:

Search with external clean function
========================================

If you want to seperate clean and search process, you can set ``clean`` in ``search_topn_matches()`` to ``False`` and use an **external clean function** :ref:`clean_spectrum() <clean-spectrum-function>` in ``ms_entropy``.

You can use the ``clean_spectrum()`` function in ``ms_entropy`` to clean the query spectrum and then use individual search functions to search the library.

.. code-block:: python

    from ms_entropy import clean_spectrum

    query_spectrum = {"precursor_mz": 150.0,
                        "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32)}

    precursor_ions_removal_da = 1.6

    query_spectrum['peaks'] = clean_spectrum(
        peaks = query_spectrum['peaks'],
        max_mz = query_spectrum['precursor_mz'] - precursor_ions_removal_da
    )

Now the query_spectrum is cleaned and ready for search. Then pass it to the ``search_topn_matches()`` with ``clean`` set to ``False``.

.. code-block:: python

    result=entropy_search.search_topn_matches(
            precursor_mz=query_spectrum['precursor_mz'],
            peaks=query_spectrum['peaks'],
            ms1_tolerance_in_da=0.01, # You can change ms1_tolerance_in_da as needed.
            ms2_tolerance_in_da=0.02, # You can change ms2_tolerance_in_da as needed.
            method='open', # Or 'neutral_loss' or 'hybrid' or 'identity'.
            precursor_ions_removal_da=1.6, # Peaks with m/z greater than ``precursor_mz - precursor_ions_removal_da`` are removed during cleaning. 
            noise_threshold=0.01, # Relative intensity threshold for noise filtering during cleaning. Peaks with intensity ``< noise_threshold * max(intensity)`` are removed.  
            min_ms2_difference_in_da=0.05, # Minimum spacing allowed between MS/MS peaks during cleaning.  
            max_peak_num=None, # Maximum number of peaks to keep after cleaning.  
            clean=False, # If you don't want to use the internal clean process in this function, set it to False.
            topn=3, # You can change topn as needed.
            need_metadata=True, # Set it to True if need metadata.
    )


You can also pass the query spectrum into the search functions mentioned in :doc:`Useful Functions </dynamic_useful_functions>` like this:

.. code-block:: python

    # Identity search
    entropy_similarity = entropy_search.identity_search(
        precursor_mz = query_spectrum['precursor_mz'],
        peaks = query_spectrum['peaks'],
        ms1_tolerance_in_da = 0.01,
        ms2_tolerance_in_da = 0.02
    )


Tools: external clean function
========================================

Both ``search_topn_matches()`` function and ``search()`` function include internal cleaning of query spectrum before performing search. 
If you want to seperate these two process, you can set ``clean`` in these two functions to ``False`` and use an external clean function. See :ref:`example <external-clean-before-search>`.

You can use the ``clean_spectrum`` function to clean the query spectrum and then use individual search functions to search the library. 

.. _clean-spectrum-function:

Clean spectrum
---------------

Before performing a spectra search, the query spectrum should be pre-processed using the ``clean_spectrum()`` function in ``ms_entropy``. This function accomplishes the following:

1. Remove empty peaks (m/z <= 0 or intensity <= 0).

2. Remove peaks with m/z values greater than ``precursor_mz - precursor_ions_removal_da`` (removes precursor ions to improve the quality of spectral comparison).

3. Centroid the spectrum by merging peaks within +/- ``min_ms2_difference_in_da`` and sort the resulting spectrum by m/z.

4. Remove peaks with intensity less than ``noise_threshold`` * maximum intensity.

5. Retain only the top max_peak_num peaks and remove all others.

6. Normalize the intensity to sum to 1.

