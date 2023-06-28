================
Useful functions
================

Reading Spectra from a File
===========================
For ease of use, a function named ``read_one_spectrum`` is provided in the ``ms_entropy`` package, allowing you to easily read spectra from a file. Here is an example of how you can use it:

.. code-block:: python

    from ms_entropy import read_one_spectrum
    for spectrum in read_one_spectrum('path/to/spectrum/file'):
        print(spectrum)

This function returns a dictionary, where each key-value pair corresponds to a specific metadata of the spectrum.

Currently, the ``read_one_spectrum`` function supports the following file formats: ``.mgf``, ``.msp``, ``.mzML``, and ``.lbm2`` from the MS-DIAL software.
