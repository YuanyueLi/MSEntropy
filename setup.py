from distutils.core import setup, Extension
from setuptools import find_packages
from Cython.Build import cythonize
import numpy as np
import os


os.environ["CFLAGS"] = "-O3 -Wno-cpp -Wno-unused-function"

setup(
    name="ms_entropy",
    version="0.6.0",
    license="Apache License 2.0",
    author="Yuanyue Li",
    url="https://github.com/YuanyueLi/SpectralEntropy",
    packages=find_packages(where=".", exclude=["tests", "docs", "examples", "manuscript", "dist", "build"]),
    python_requires=">=3.8",
    install_requires=[
        "numpy >= 1.18",
        "cython >= 0.29",
    ],
    extras_require={"gpu": ["cupy >= 12.0.0"]},
    keywords=[
        "ms entropy",
        "ms spectral entropy",
        "spectral entropy",
        "spectral similarity",
        "entropy similarity",
        "entropy",
        "entropy search",
        "flash entropy",
        "flash entropy similarity",
        "flash entropy search",
    ],
    package_dir={"": "."},
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
