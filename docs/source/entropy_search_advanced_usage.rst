==============
Advanced usage
==============


Run Flash entropy search with limited memory
============================================

This method is useful when you are dealing with a very large spectral library and your computer's memory is limited.

To achieve this, while constructing the ``FlashEntropySearch`` object, you need to set the ``path_data`` parameter to the path of the index file, and set the ``low_memory`` parameter to ``True``. Then read the pre-built index file by calling the ``read`` method. After that, the rest of the code is the same as usual.

.. code-block:: python

    from ms_entropy import FlashEntropySearch

    # Instead of using this:
    # entropy_search = FlashEntropySearch()
    # Use this:
    entropy_search = FlashEntropySearch(path_data='path/to/library/index', low_memory=True)
    entropy_search.read()

    # Then the reset of the code is the same as usual.
    # entropy_search.search(...)
    # ...... (the reset of the code is the same as usual)

The index built in normal mode and low memory mode is identical. If you use our ``write`` and ``read`` methods to save and load the index, you can use the index in normal mode and low memory mode interchangeably. For example, you can build the index in normal mode, save it to disk with the ``write`` method. After that, you can initialize the ``FlashEntropySearch`` object with ``path_data`` parameter which points to the index file, and set ``low_memory`` parameter to ``True``, then call the ``read`` method to load the index, and proceed with the search as usual.


Run Flash entropy search with multiple cores
============================================

When you have many query spectra and your computer has multiple cores, you can use these multiple cores to speed up the search. Python's built-in ``multiprocessing`` module can be utilized for this purpose.

To avoid the overhead of initializing the ``FlashEntropySearch`` object in each process, you can use the ``save_memory_for_multiprocessing`` method. This function copies the index data to shared memory, avoiding the overhead of copying the index data to each process. You can then use the ``initializer`` and ``initargs`` parameters of the ``Pool`` class to initialize the ``FlashEntropySearch`` object in each process.

Here's an example on how to use the ``multiprocessing`` module to calculate the maximum entropy similarity of 100 query spectra with 4 cores.

.. code-block:: python
    
    import multiprocessing as mp
    
    query_spectrum = {"precursor_mz": 150.0,
                      "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32)}

    # Let's say you have 100 query spectra.
    query_spectra_list = [query_spectrum] * 100

    # And you want to use 4 cores to speed up the search.
    THREADS = 4

    entropy_search = FlashEntropySearch()
    entropy_search.build_index(spectral_library)

    # Call this function to save memory for multiprocessing.
    entropy_search.save_memory_for_multiprocessing()

    def func_max_similarity(precursor_mz, peaks):
        entropy_search = func_max_similarity.entropy_search
        similarity = entropy_search.search(precursor_mz=precursor_mz, peaks=peaks, method="neutral_loss")
        return np.max(similarity["neutral_loss_search"])

    def init_worker(entropy_search, ):
        func_max_similarity.entropy_search = entropy_search
        return None

    pool = mp.Pool(THREADS, initializer=init_worker, initargs=(entropy_search, ))
    max_entropy_similarity = pool.starmap(func_max_similarity, [(spectrum["precursor_mz"], spectrum["peaks"]) for spectrum in query_spectra_list])

.. note:: 
    When using the multiple cores, you always need to keep in mind of the memory usage.
    
    For instance, if you're attempting to search 1,000,000 spectra against a spectral library containing 1,000,000 spectra, the result similarity matrix will be **4*1,000,000*1,000,000 = 3.6 TB!** Therefore, returning the entire similarity matrix may not be practical.
    
    Instead, you can process the result similarity and only return the processed result. This saves a lot of memory and computational time needed for copying large amounts of memory. For example, if you only return the top 10 similarities for each query spectrum, the memory usage will be **4*1,000,000*10 = 38 MB**, which is significantly more efficient.

Run Flash entropy search on GPU
===============================

When you have a GPU and searching a single spectrum takes more than 0.1 seconds, you can use the GPU to speed up the search. To do this, you'll need to install the `Cupy <https://cupy.dev/>`_ package first. You can then use the ``target`` parameter set to ``gpu`` to use the GPU.

.. code-block:: python

    from ms_entropy import FlashEntropySearch
    entropy = FlashEntropySearch()
    entropy_search.build_index(spectral_library)

    # Instead of using this:
    # entropy_similarity = entropy_search.search(
    #     precursor_mz = 150.0,
    #     peaks = [[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]],
    # )

    # Use this:
    entropy_similarity = entropy_search.search(
        precursor_mz = 150.0,
        peaks = [[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]],
        target = 'gpu'
    )

    # The rest of your code remains the same.

The return values when calculating with ``CPU`` and ``GPU`` are the same. Hence, you can use the same code to process the result. Running computations on a GPU can substantially speed up your program if you're performing large-scale spectral library searching.
