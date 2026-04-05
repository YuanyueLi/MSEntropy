.. MS Entropy documentation master file, created by
   sphinx-quickstart on Sat Jun  3 23:52:16 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction to MS Entropy
==========================

Thank you for using the MS Entropy package, you can download the code from `MS Entropy Github repository <https://github.com/YuanyueLi/MSEntropy>`_.
If you encounter any issues, queries or need support, don't hesitate to contact :email:`Yuanyue Li <mail@yli.one>`

If you find this package useful, please consider citing the following papers:

Li, Y., Fiehn, O., Flash entropy search to query all mass spectral libraries in real time.  *Nat Methods* **20**, 1475-1478 (2023). `https://doi.org/10.1038/s41592-023-02012-9 <https://doi.org/10.1038/s41592-023-02012-9>`_

Li, Y., Kind, T., Folz, J. _et al._ Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. *Nat Methods* **18**, 1524-1531 (2021). `https://doi.org/10.1038/s41592-021-01331-z <https://doi.org/10.1038/s41592-021-01331-z>`_


This package contains the following modules:

1.  **Classical entropy functions:** These functions are ideal for calculating the classical spectral entropy and entropy similarity for a small set of spectrum pairs.


2.  **Flash Entropy Search Algorithm Functions:** If your task involves searching spectra against a large spectral library, this algorithm will significantly improve the speed of the search process.

3.  **Dynamic Entropy Search Algorithm Functions:** This module supports fast index build, incremental updates of reference spectra library, and spectra search against a large spectral library with high efficiency.

4.  **Repository Search Functions:** This module supports search and personalized incremental updates of our prebuilt index. This prebuilt resource comprises all publicly available spectra (1,473,028,876 MS/MS spectra from MassIVE/GNPS, MetaboLights, and Metabolomics Workbench), enabling global spectra comparisons.

.. toctree::
   :maxdepth: 1

   self
   install

.. toctree::
   :maxdepth: 1
   :caption: 📊 Classical entropy functions 📊

   classical_spectral_entropy
   classical_entropy_similarity
   classical_useful_functions
   classical_api

.. toctree::
   :maxdepth: 1
   :caption: ⚡ Flash Entropy Search ⚡

   entropy_search_quickstart
   entropy_search_basic_usage
   entropy_search_advanced_usage
   entropy_search_useful_functions
   entropy_search_examples
   entropy_search_api

.. toctree::
   :maxdepth: 1
   :caption: ⚡ Dynamic Entropy Search ⚡

   dynamic_quickstart
   dynamic_basic_usage_index_construction
   dynamic_basic_usage_search
   dynamic_useful_functions
   dynamic_api

.. toctree::
   :maxdepth: 1
   :caption: 🌐 Repository Search 🌐

   repository_search_basic_usage
   repository_api