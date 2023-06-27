============
Installation
============

To use this package, you will need to have ``Python >= 3.7`` installed on your system.


From PyPI
============

To install the latest version of the package from PyPI, run the following command:

.. code-block:: bash

  pip install ms_entropy


Before installing the package, you will need to have ``C compiler and Python development headers.`` installed on your system:

  - On Linux, you may need to install the ``gcc`` and ``python-dev`` (for apt) or ``python-devel`` (for yum) packages first.
  - On Windows, you will need to install the `Microsoft Visual C++ 14.0 or greater Build Tools <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_ first.
  - On macOS, you will need to install the ``Xcode Command Line Tools`` first.



From source
============

This package can run in pure Python mode, or in Cython mode (which will compile the Cython code to run faster).

Pure python mode
----------------

When run in pure Python mode, you only need copy the folder ``ms_entropy`` to your project directory. Or set your ``PYTHONPATH`` environment variable to include the ``ms_entropy`` folder.


Cython mode
-----------

In Cython mode, you will need to compile the Cython code first. To do this, clone the repository and run the following commands:

Requirements
^^^^^^^^^^^^

Before compiling, you will need to have ``C compiler and Python development headers.`` installed on your system:

  - On Linux, you may need to install the ``gcc`` and ``python-dev`` (for apt) or ``python-devel`` (for yum) packages first.
  - On Windows, you will need to install the `Microsoft Visual C++ 14.0 or greater Build Tools <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_ first.
  - On macOS, you will need to install the ``Xcode Command Line Tools`` first.

Downlaod and compile
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  git clone https://github.com/YuanyueLi/MSEntropy.git
  cd MSEntropy
  python setup.py build_ext --inplace
  
Then you can copy the ``ms_entropy`` folder to your project directory. Or set your ``PYTHONPATH`` environment variable to include the ``ms_entropy`` folder.


Dependencies
^^^^^^^^^^^^

The package depends on the following Python packages, which will be installed automatically when you install the package from PyPI. If you install the package from source, you will need to install them manually:


- ``Numpy >= 1.9.3``

- ``Cython >= 0.26.1`` (optional, which will compile the Cython code to run faster)

- ``Cupy >= 12.0.0`` (optional, required for GPU acceleration)

- ``pyteomics >= 4.6`` (optional, required for reading mzML files)

- ``lz4 >= 4.3.2``, ``msgpack >= 1.0.5`` (optional, required for reading .lbm2 files format from MS-DIAL)


The code is tested on Windows, Linux and macOS.


Test
====

To test that the package is working correctly, run the example.py script:


.. code-block:: bash

  cd examples
  python example.py
