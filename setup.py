from setuptools import find_packages, setup, Extension
from distutils.util import convert_path
from Cython.Build import cythonize
import numpy as np
import os


os.environ["CFLAGS"] = "-O3 -Wno-cpp -Wno-unused-function"

# Read the version from the package file
main_ns = {}
ver_path = convert_path("ms_entropy/version.py")
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)


setup(
    name="ms_entropy",
    version=main_ns["__version__"],
    url="https://github.com/YuanyueLi/MSEntropy",
    license="Apache License 2.0",
    license_file="LICENSE",
    platforms=["Linux", "MacOS", "Windows"],
    packages=find_packages(where='.', exclude=['tests', 'docs', 'examples', 'manuscript', 'dist', 'build']),
    package_dir={"": "."},

    python_requires=">=3.7",
    install_requires=[
        "numpy >= 1.9.3",
    ],
    extras_require={
        "all": ["lz4 >= 4.3.2", "msgpack >= 1.0.5", "pyteomics >= 4.6"],
        "gpu": ["cupy >= 12.0.0"]
    },

    ext_modules=cythonize(
        [
            Extension(
                "ms_entropy.spectra.entropy_cython",
                [r"ms_entropy/spectra/entropy_cython.pyx", r"ms_entropy/spectra/CleanSpectrum.c", r"ms_entropy/spectra/SpectralEntropy.c"],
            ),
            Extension(
                "ms_entropy.entropy_search.fast_flash_entropy_search_cpython",
                [r"ms_entropy/entropy_search/fast_flash_entropy_search_cpython.pyx"],
            ),
        ],
        annotate=False,
        compiler_directives={
            "language_level": "3",
            "cdivision": True,
            "boundscheck": False,  # turn off bounds-checking for entire function
            "wraparound": False,  # turn off negative index wrapping for entire function
        },
    ),
    include_dirs=[np.get_include()],
)

# python setup.py sdist
# python setup.py build_ext --inplace
