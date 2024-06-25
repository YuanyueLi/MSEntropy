# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

sys.path.insert(0, os.path.abspath("../../"))
sys.path.insert(0, os.path.abspath("../../ms_entropy"))

from ms_entropy import __version__

project = "MS Entropy"
copyright = "2023, Yuanyue Li"
author = "Yuanyue Li"
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", 
              "sphinx.ext.doctest", 
              "sphinx.ext.autosummary", 
              "sphinx.ext.viewcode", 
              "sphinx.ext.githubpages",
              "sphinxcontrib.email",
              "numpydoc"]

# apidoc_module_dir = "../../ms_entropy"
# apidoc_output_dir = "./api"
# apidoc_separate_modules = True
# apidoc_module_first = True
numpydoc_show_class_members = False
autoclass_content = 'both'

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
