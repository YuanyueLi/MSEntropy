============
Installation
============

MS Entropy requires ``Python >= 3.7`` installed on your system. It has been tested on Windows, Linux, and macOS platforms.


Dependencies
============

If you install MS Entropy from PyPI, the necessary dependencies will be installed automatically. However, if you are installing from the source, install these manually:

- ``Numpy >= 1.9.3``
- ``Cython >= 0.26.1`` (Optional: required for Cython mode)
- ``Cupy >= 12.0.0`` (Optional: required for GPU acceleration)
- ``pyteomics >= 4.6`` (Optional: required for reading mzML files)
- ``lz4 >= 4.3.2``, ``msgpack >= 1.0.5`` (Optional: required for reading .lbm2 file format from MS-DIAL)


Installing from PyPI
====================

Before starting the installation process, ensure that you have a ``C compiler and Python development headers`` installed:

- On Linux, install the ``gcc`` and ``python-dev`` (for apt) or ``python-devel`` (for yum) packages.
- On Windows, install the `Microsoft Visual C++ 14.0 or greater Build Tools <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_.
- On macOS, install the ``Xcode Command Line Tools``.

To install the latest version of MS Entropy from PyPI, use the following command:

.. code-block:: bash

  pip install ms_entropy


If you need to use **GPU** acceleration, you need to install the latest version of CuPy first (recommended). Or you can install the latest version of MS Entropy with GPU support from PyPI:

.. code-block:: bash

  pip install ms_entropy[gpu]


Installing from Source
======================

MS Entropy can operate in pure Python mode or in Cython mode for a performance boost.

Pure Python Mode
----------------

In this mode, simply copy the ``ms_entropy`` folder to your project directory or include the ``ms_entropy`` directory in your ``PYTHONPATH`` environment variable.

Cython mode
-----------

To operate in Cython mode, compile the Cython code after cloning the repository:

- Requirements

  Before compiling, you will need to have ``C compiler and Python development headers.`` installed on your system:

  - On Linux, you may need to install the ``gcc`` and ``python-dev`` (for apt) or ``python-devel`` (for yum) packages first.
  - On Windows, you will need to install the `Microsoft Visual C++ 14.0 or greater Build Tools <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_ first.
  - On macOS, you will need to install the ``Xcode Command Line Tools`` first.

- Downlaod and compile

  .. code-block:: bash

    git clone https://github.com/YuanyueLi/MSEntropy.git
    cd MSEntropy
    python setup.py build_ext --inplace
    
  Then, as in Pure Python Mode, you can either copy the ms_entropy folder to your project directory or include it in your PYTHONPATH environment variable.


Testing
=======

To verify the successful installation and operation of the package, run the example.py script:


.. code-block:: bash

  cd examples
  python example.py

This script will help confirm that the package is working as expected.
