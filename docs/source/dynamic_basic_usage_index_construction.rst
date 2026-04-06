=================================
Basic usage - Index Construction
=================================

In order to perform spectra comparison more efficiently, spectral data should be loaded into a reference library in the form of index.


Step 0: Prepare the library spectra
===================================

Suppose you have a lot of spectra and want to build library based on them, you need to format them like this:

.. code-block:: python

    import numpy as np
    # For each spectral library, it is a list consisting of multiple dictionaries of MS2 spectra.

    # For each spectrum, 'precursor_mz' and 'peaks' are necessary. 
    # 'precursor_mz' should be a float, and 'peaks' should be a 2D np.ndarray like np.ndarray([[m/z, intensity], [m/z, intensity], [m/z, intensity]...], dtype=np.float32).


    spectra_1_for_library = [{
        "id": "Demo spectrum 1",
        "precursor_mz": 150.0,
        "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [103.0, 1.0]], dtype=np.float32), 
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
        "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [302.0, 1.0]], dtype=np.float32),
    },
        ]

    spectra_2_for_library ... # Similar to spectra_1_for_library
    spectra_3_for_library ... # Similar to spectra_1_for_library


The keys ``precursor_mz`` and ``peaks`` are necessary for this format.
Other keys are optional and are considered as metadata of the spectrum, helping the identification.

Then you can have your spectra lists to be added into the library.


Step 1: Perform updates
=======================

Initial construction
----------------------
Suppose that you want to construct an index with spectra_1_for_library at first:


.. code-block:: python

    # Firstly, import DynamicEntropySearch.
    from ms_entropy import DynamicEntropySearch

    # Secondly, assign the path for your library.
    entropy_search=DynamicEntropySearch(
            path_data=path_of_your_library, 
            max_ms2_tolerance_in_da=0.024, # Maximum MS/MS tolerance (in Daltons) used during spectrum search.
            extend_fold=3, # Expansion factor for preallocated storage in each m/z block. Determines ``reserved_len = data_len * extend_fold``. 
            mass_per_block=0.05, # m/z step size for creating the index blocks.
            num_per_group=100_000_000, # Number of spectra assigned to each group. 
            cache_list_threshold=1_000_000, # Number of spectra to accumulate in memory before writing them to disk.
            max_indexed_mz=1500.00005, # Maximum m/z value to index. Ions above this threshold are grouped into a single block. 
            intensity_weight="entropy",  # "entropy" or None.Determines whether intensities are entropy-weighted. 
    )

    # Thirdly, add spectra list into the library one by one. By default, there will be a built-in cleaning function.
    entropy_search.add_new_spectra(spectra_list=spectra_1_for_library)

    # Lastly, call build_index() and write() to end the adding operation.
    entropy_search.build_index()
    entropy_search.write()


There are some tips in this process:

- It is necessary to initialize ``DynamicEntropySearch`` using a specified ``path_data``, which is the path of your library. The reset of the parameters are optional. If it is a new library, the value of ``path_data`` should be new and will be created in the initialization of class.


- If you only want to build index for open search, you can set ``index_for_neutral_loss`` in ``add_new_spectra()`` and ``build_index()`` to ``False``. However, after doing this, you couldn't perform neutral loss search or hybrid search under this library anymore. Besides, adding neutral loss mass index into this library is violated too. This means that once the ``index_for_neutral_loss`` in the ``add_new_spectra()`` function as well as ``build_index()`` function are set to ``False``, they should remain ``False`` from then on. Any violation can cause errors.


- There is a built-in cleaning function in ``add_new_spectra()``. Peaks can be cleaned using this function. Cleaning is a necessary part in construction. If ``clean`` in ``add_new_spectra()`` is set to ``False``, use external clean function :ref:`clean_spectrum() <clean-spectrum-function>` in ``ms_entropy`` to process spectra. See the following example. Lack of cleaning can lead to error.


- It is necessary to call ``build_index()`` and ``write()`` lastly after all ``add_new_spectra()`` as the end of adding operation to make sure all the spectra are loaded into the index.

.. _external-clean-before-add:

If use **external clean function** to process spectra in construction:

.. code-block:: python

    # Firstly, import.
    from ms_entropy import DynamicEntropySearch
    from ms_entropy import clean_spectrum

    # Secondly, assign the path for your library.
    entropy_search=DynamicEntropySearch(
            path_data=path_of_your_library, 
            max_ms2_tolerance_in_da=0.024, # Maximum MS/MS tolerance (in Daltons) used during spectrum search.
            extend_fold=3, # Expansion factor for preallocated storage in each m/z block. Determines ``reserved_len = data_len * extend_fold``. 
            mass_per_block=0.05, # m/z step size for creating the index blocks.
            num_per_group=100_000_000, # Number of spectra assigned to each group. 
            cache_list_threshold=1_000_000, # Number of spectra to accumulate in memory before writing them to disk.
            max_indexed_mz=1500.00005, # Maximum m/z value to index. Ions above this threshold are grouped into a single block. 
            intensity_weight="entropy",  # "entropy" or None.Determines whether intensities are entropy-weighted. 
    )

    # Manually clean using external clean function
    precursor_ions_removal_da = 1.6

    spectra_1_for_library_clean=[]
    for spec in spectra_1_for_library:
        spec['peaks'] = clean_spectrum(
            peaks = spec['peaks'],
            max_mz = spec['precursor_mz'] - precursor_ions_removal_da)
        if len(spec['peaks']) > 0:
            spectra_1_for_library_clean.append(spec)

    # Thirdly, add spectra list into the library one by one. Set `clean` to `False` because spectra have been cleaned before.
    entropy_search.add_new_spectra(spectra_list=spectra_1_for_library_clean, clean=False)

    # Lastly, call build_index() and write() to end the adding operation.
    entropy_search.build_index()
    entropy_search.write()


Generally, we **recommend internal clean** in ``add_new_spectra()``.

Note that three parameters:

(1) ``max_ms2_tolerance_in_da`` in the initialization of class ``DynamicEntropySearch()``

(2) ``min_ms2_difference_in_da`` in ``add_new_spectra()``

(3) ``ms2_tolerance_in_da`` in search functions of ``DynamicEntropySearch()``

should follow this rule: ``min_ms2_difference_in_da`` > ``max_ms2_tolerance_in_da`` * 2 >= ``ms2_tolerance_in_da`` * 2.

An error will be reported if the condition is not met.

Once these steps are complete, you will find a folder, which serves as the library, at the ``path_data``. In this folder, several binary files and one or more subfolders can be found. These binary files record the information of subfolders and metadata. Each subfolder refers to a group — the organizational unit directly under a library. These subfolders are numerically named, starting from 0.

Example structure — one library containing 3 groups:

.. code-block:: bash

    path_of_your_library/
    ├── 0/
    ├── 1/
    ├── 2/
    ├── group_start.pkl
    ├── metadata_start_loc.bin
    └── metadata.pkl 


The library with index of spectra_1_for_library is created. The spectra_1_for_library is now saved as index in this library in group 0. You can fetch the library whenever you want.


Subsequent construction
------------------------

Next time, if you want to update the index with spectra_2_for_library and spectra_3_for_library, just select the correct path and execute the update:

.. code-block:: python

    # Firstly, import DynamicEntropySearch.
    from ms_entropy import DynamicEntropySearch

    # Secondly, choose the existing library with corresponding path. This library is built with spectra_1_for_library last time.
    entropy_search=DynamicEntropySearch(
            path_data=path_of_your_library, 
            max_ms2_tolerance_in_da=0.024, # Maximum MS/MS tolerance (in Daltons) used during spectrum search.
            extend_fold=3, # Expansion factor for preallocated storage in each m/z block. Determines ``reserved_len = data_len * extend_fold``. 
            mass_per_block=0.05, # m/z step size for creating the index blocks.
            num_per_group=100_000_000, # Number of spectra assigned to each group. 
            cache_list_threshold=1_000_000, # Number of spectra to accumulate in memory before writing them to disk.
            max_indexed_mz=1500.00005, # Maximum m/z value to index. Ions above this threshold are grouped into a single block. 
            intensity_weight="entropy",  # "entropy" or None.Determines whether intensities are entropy-weighted. 
    )

    # Thirdly, add spectra list into the library one by one.
    entropy_search.add_new_spectra(spectra_list=spectra_2_for_library)
    entropy_search.add_new_spectra(spectra_list=spectra_3_for_library)

    # Lastly, call build_index() and write() to end the adding operation.
    entropy_search.build_index()
    entropy_search.write()


Now the library in ``path_of_your_library`` contains index of 3 spectra lists: spectra_1_for_library, spectra_2_for_library and spectra_3_for_library. They may be distributed in one or more groups, depending on the number of spectra and the value of ``num_per_group``.


Tools: external clean function
========================================

Function ``add_new_spectra()`` includes internal cleaning of reference spectra before actually constructing index. 

If you want to seperate these two process, you can set ``clean`` in this function to ``False`` and use an external clean function. See :ref:`example <external-clean-before-add>`.

You can use the ``clean_spectrum()`` function to clean the reference spectra and then use ``add_new_spectra()`` to perform updates. 

.. _clean-spectrum-function:

Clean spectrum
---------------

Before performing a library update, the reference spectrum should be pre-processed using the ``clean_spectrum()`` function in ``ms_entropy``. This function accomplishes the following:

1. Remove empty peaks (m/z <= 0 or intensity <= 0).

2. Remove peaks with m/z values greater than ``precursor_mz - precursor_ions_removal_da`` (removes precursor ions to improve the quality of spectral comparison).

3. Centroid the spectrum by merging peaks within +/- ``min_ms2_difference_in_da`` and sort the resulting spectrum by m/z.

4. Remove peaks with intensity less than ``noise_threshold`` * maximum intensity.

5. Retain only the top max_peak_num peaks and remove all others.

6. Normalize the intensity to sum to 1.

