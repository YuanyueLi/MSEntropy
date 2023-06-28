.. MS Entropy documentation master file, created by
   sphinx-quickstart on Sat Jun  3 23:52:16 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction to MS Entropy
==========================

Thank you for using the MS Entropy package, you can download the code from `MS Entropy Github repository <https://github.com/YuanyueLi/MSEntropy>`_.
If you encounter any issues, queries or need support, don't hesitate to contact :email:`Yuanyue Li <liyuanyue@gmail.com>`

This package comprises the following modules:

1.  **Classical entropy functions:** These functions are ideal for calculating the classical spectral entropy and entropy similarity for a small set of spectrum pairs.

3.  **Flash Entropy Search Algorithm Functions:** If your task involves searching spectra against a large spectral library, this algorithm will significantly improve the speed of the search process.


For a comprehensive understanding of the package, refer to our documentation hosted on `Read the Docs <https://msentropy.readthedocs.io>`_.
Additionally, you may find practical examples in the `examples <https://github.com/YuanyueLi/MSEntropy/tree/main/examples>`_ folder on our Github repository.


.. toctree::
   :maxdepth: 2

   self
   install
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: 📊 Classical entropy functions:

   classical_useful_functions
   classical_api

.. toctree::
   :maxdepth: 2
   :caption: ⚡ Flash Entropy Search:

   entropy_search_basic_usage
   entropy_search_advanced_usage
   entropy_search_useful_functions
   entropy_search_examples
   entropy_search_api

