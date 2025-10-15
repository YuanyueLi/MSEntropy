============
Installation
============

**MS Entropy** requires **Python >= 3.8**. While older versions like Python 3.7 *may* work, they are not guaranteed to be compatible.  
The package has been tested on **Windows**, **Linux**, and **macOS**.


Installing from PyPI
====================

To install the latest version of MS Entropy from PyPI, use the following command:

.. code-block:: bash

  pip install ms_entropy

MS Entropy includes many useful features, such as support for reading ``.mzML``, ``.msp``, ``.mgf``, and ``.lbm2`` files.  
These features are optional and not included by default to keep the package lightweight.  
To install MS Entropy with all optional features, use:


.. code-block:: bash

  pip install ms_entropy[all]
GPU Support
-----------

To enable **GPU acceleration**, it's recommended to install the latest version of **CuPy** separately.  
Alternatively, you can install MS Entropy with GPU support via:

.. code-block:: bash

  pip install ms_entropy[gpu]

Full Functionality (GPU + Extras)
---------------------------------

To install MS Entropy with both GPU support and all optional features:

.. code-block:: bash

  pip install ms_entropy[gpu,all]


Dependencies
============

If you install MS Entropy from PyPI, the necessary dependencies will be installed automatically. However, if you are installing from the source, install these manually:

- ``Numpy >= 1.9.3``
- ``Cython >= 0.26.1`` (Optional: required for Cython mode)
- ``Cupy >= 12.0.0`` (Optional: required for GPU acceleration)
- ``pyteomics >= 4.6`` (Optional: required for reading mzML files)
- ``lz4 >= 4.3.2``, ``msgpack >= 1.0.5`` (Optional: required for reading .lbm2 file format from MS-DIAL)


Installing from Source
======================

MS Entropy can operate in either **pure Python mode** or **Cython mode** for performance improvements.

Pure Python Mode
----------------

In this mode, simply copy the ``ms_entropy`` folder to your project directory,  
or include it in your ``PYTHONPATH`` environment variable.

Cython mode
-----------

To operate in Cython mode, compile the Cython code after cloning the repository:

Requirements
~~~~~~~~~~~~

  Before compiling, you will need to have ``C compiler and Python development headers.`` installed on your system:

  - On Linux, you may need to install the ``gcc`` and ``python-dev`` (for apt) or ``python-devel`` (for yum) packages first.
  - On Windows, you will need to install the `Microsoft Visual C++ 14.0 or greater Build Tools <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_ first.
  - On macOS, you will need to install the ``Xcode Command Line Tools`` first.

Download and Compile
~~~~~~~~~~~~~~~~~~~~

  .. code-block:: bash

    git clone https://github.com/YuanyueLi/MSEntropy.git
    cd MSEntropy
    python setup.py build_ext --inplace
    
After compilation, you can either copy the ``ms_entropy`` folder to your project directory  
or include it in your ``PYTHONPATH``.


Testing
=======

To verify the successful installation and operation of the package, run the example.py script:


.. code-block:: bash

  cd examples
  python example.py

This script will help confirm that the package is working as expected.
