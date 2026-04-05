===========
Basic usage
===========


Overview
========

We offer prebuilt index for metabolomics repositories, comprising more than 1.4 billion spectra.
Users can use ``RepositorySearch`` class to search against these public metabolomics repositories. 
We have built the indexes and uploaded them to `Hugging Face repository <https://huggingface.co/datasets/YuanyueLiZJU/dynamic_entropy_search/tree/main>`_.
By downloading these indexes and extracting them to a specified path on your local machine, you can perform search like this:

Firstly, assign the path of the prebuilt indexes as the ``path_data`` of ``RepositorySearch`` class.

Remember to prepare query spectrum in correct format and clean it (see aforementioned points to prepare the format in :doc:`Dynamic Entropy Search: Basic usage - Search </dynamic_basic_usage_search>`). Note that key ``charge`` here is **necessary** too. Set it to 1 or -1.

Then perform search, and you can get top few results.

.. code-block:: python

    from ms_entropy import RepositorySearch
    import numpy as np
    from ms_entropy import clean_spectrum

    # Instantiation
    entropy_search=RepositorySearch(path_data=path_repository_indexes)

    # Prepare query spectra
    precursor_ions_removal_da=1.6

    query_spec_1={
        "peaks":np.array([[217.0, 1.5], [234.0, 0.8], [398.0, 2.0]]),
        "precursor_mz":455.0,
        "charge":-1                
    }

    query_spec_2={
        "peaks":np.array([[123.0, 1.0], [126.0, 0.7], [101.0, 4.0]]),
        "precursor_mz":250.0,
        "charge":1                
    }

    query_spec_3={
        "peaks":np.array([[200.0, 0.9], [101.0, 3.2], [202.0, 1.7]]),
        "precursor_mz":345.0,
        "charge":-1                
    }
    query_spectra=[query_spec_1, query_spec_2, query_spec_3]

    # Clean query spectra
    for query_spec in query_spectra:
        query_spec['peaks']=clean_spectrum(
                peaks=query_spec['peaks'],
                max_mz = query_spec['precursor_mz'] - precursor_ions_removal_da
            )
        
    # Perform search and output results
    for i, query_spec in enumerate(query_spectra):
        result=entropy_search.search_topn_matches(
            method="open", 
            charge=query_spec['charge'],
            precursor_mz=query_spec['precursor_mz'],
            peaks=query_spec['peaks'],
            topn=3 # can be changed
            )

        print(f"Query spectrum {i} matches:{result}\n")

An example result:


.. code-block:: bash

    Query spectrum 0 matches:[{'file_name': 'gnps/MSV000080555/A9_RA9_01_8358.mzML.gz', 'scan': np.uint64(1074), 'similarity': np.float64(0.758648693561554), 'spec_idx': np.uint64(187632338)}, {'file_name': 'gnps/MSV000094528/20230131_pluskal_mce_1D2_H4_id_negative.mzML.gz', 'scan': np.uint64(319), 'similarity': np.float64(0.6575593948364258), 'spec_idx': np.uint64(231674808)}, {'file_name': 'metabolights/MTBLS700/ns154.mzML.gz', 'scan': np.uint64(1664), 'similarity': np.float64(0.6054909229278564), 'spec_idx': np.uint64(261307506)}]

    Query spectrum 1 matches:[{'file_name': 'gnps/MSV000079098/Stairs_Gz140_57H3_RH3_01_651.mzML.gz', 'scan': np.uint64(780), 'similarity': np.float64(0.8724607229232788), 'spec_idx': np.uint64(168727792)}, {'file_name': 'gnps/MSV000079098/Roff_SEA12_64F2_GF2_01_721.mzML.gz', 'scan': np.uint64(2860), 'similarity': np.float64(0.8721256852149963), 'spec_idx': np.uint64(249982900)}, {'file_name': 'gnps/MSV000080141/A2154D_34E05_RE5_01_25781.mzML.gz', 'scan': np.uint64(487), 'similarity': np.float64(0.8716126680374146), 'spec_idx': np.uint64(171140687)}]

    Query spectrum 2 matches:[{'file_name': 'metabolomics_workbench/ST003745/x01997_NEG.mzML.gz', 'scan': np.uint64(1139), 'similarity': np.float64(0.7167922258377075), 'spec_idx': np.uint64(259794461)}, {'file_name': 'gnps/MSV000090773/B2_EC_P_NEG_QC_Cond_MSMS_AutoMSMS_1.mzML.gz', 'scan': np.uint64(421), 'similarity': np.float64(0.6405000686645508), 'spec_idx': np.uint64(325131566)}, {'file_name': 'gnps/MSV000090773/B2_EC_P_NEG_QC_Cond_MSMS_AutoMSMS_1.mzML.gz', 'scan': np.uint64(423), 'similarity': np.float64(0.6405000686645508), 'spec_idx': np.uint64(325131568)}]

In this result:

- Every query spectrum gets 3 matched reference spectra because ``topn`` in ``search_topn_matches()`` is set to ``3``. This value can be changed based on your need.
- Every result contains ``file_name``, ``similarity``, ``scan`` that can be used in further identification.

By using ``RepositorySearch``, you can easily search spectra against global spectra library.


.. note::

    It should be noted that, query spectrum in this part should have a key ``charge`` with an integer value (-1 or 1). 


------------


Example
========
Another example can be found in the root directory of DynamicEntropySearch `Github <https://github.com/2bereal-me/DynamicEntropySearch>`_.

example_repository.py
------------------------

This script demonstrates how to use the RepositorySearch module.
If you want to use our prebuilt indexes, replace the ``path_data`` with the index path that you have downloaded and extracted to your local machine. Then search can be performed directly.

If you want to build or add your own reference spectra, you can follows the instruction in this example.

