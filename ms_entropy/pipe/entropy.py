import inspect
import re

import numpy as np
import pandas as pd

from ..spectra import (
    apply_weight_to_intensity,
    calculate_entropy_similarity,
    calculate_spectral_entropy,
    calculate_unweighted_entropy_similarity,
    clean_spectrum,
)
from .pipe import assign_and_apply, pipe_function

__all__ = ["add_spectral_entropy", "add_entropy_quality", "add_entropy_similarity"]


@pipe_function
def add_entropy_similarity(df, col_peaks=["peaks_a", "peaks_b"], col_precursor_mz="precursor_mz", col_result="entropy_similarity", **kwargs):
    """
    Add the entropy similarity column to the dataframe. The entropy similarity is calculated as the similarity between the two spectra in the ``col_peaks`` column.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to be added.

    col_peaks : list
        The column names of the peaks columns. Should be a list of length 2.

    col_precursor_mz : str or list or None
        The column name of the precursor m/z column.

    col_result : str
        The column name of the entropy similarity column.

    kwargs : dict
        The keyword arguments for `calculate_entropy_similarity`.

        _

    Returns
    -------
    df : pandas.DataFrame
        The dataframe with the entropy similarity column added.
    """
    precursor_mz_removal = 1.6
    if col_precursor_mz is None:
        func_search = lambda x: calculate_entropy_similarity(x[col_peaks[0]], x[col_peaks[1]], **kwargs)
    elif isinstance(col_precursor_mz, str):
        func_search = lambda x: calculate_entropy_similarity(x[col_peaks[0]], x[col_peaks[1]], max_mz=x[col_precursor_mz] - precursor_mz_removal, **kwargs)
    elif len(col_precursor_mz) == 1:
        func_search = lambda x: calculate_entropy_similarity(x[col_peaks[0]], x[col_peaks[1]], max_mz=x[col_precursor_mz[0]] - precursor_mz_removal, **kwargs)
    else:
        func_search = None

    if func_search is not None:
        similarity = df.apply(func_search, axis=1)
    else:
        peaks_a_cleaned = df.apply(lambda x: clean_spectrum(x[col_peaks[0]], max_mz=x[col_precursor_mz[0]] - precursor_mz_removal, **kwargs), axis=1)
        peaks_b_cleaned = df.apply(lambda x: clean_spectrum(x[col_peaks[1]], max_mz=x[col_precursor_mz[1]] - precursor_mz_removal, **kwargs), axis=1)
        similarity = pd.DataFrame({"peaks_a": peaks_a_cleaned, "peaks_b": peaks_b_cleaned}).apply(
            lambda x: calculate_entropy_similarity(x["peaks_a"], x["peaks_b"], **kwargs, clean_spectra=False), axis=1
        )

    return df.assign(**{col_result: similarity})


@pipe_function
def add_entropy_quality(df, col_peaks="peaks", col_result="quality", **kwargs):
    """
    Add the entropy quality column to the dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to be added.

    col_peaks : str
        The column name of the peaks column.

    col_result : str
        The column name of the quality column.

    kwargs : dict
        The keyword arguments for `calculate_spectral_entropy`.

        _

    Returns
    -------
    df : pandas.DataFrame
        The dataframe with the entropy quality column added.
    """
    entropy = df.apply(lambda x: calculate_spectral_entropy(x[col_peaks], **kwargs), axis=1)
    normalized_entropy = entropy / np.log(df[col_peaks].apply(lambda x: len(x)))
    normalized_entropy[entropy == 0] = 0
    normalized_entropy = 1 - normalized_entropy**4
    return df.assign(**{col_result: normalized_entropy})


@pipe_function
def add_spectral_entropy(df, col_peaks="peaks", col_result="entropy", **kwargs):
    """
    Add the spectral entropy column to the dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to be added.

    col_peaks : str
        The column name of the peaks column.

    col_result : str
        The column name of the entropy column.

    kwargs : dict
        The keyword arguments for `calculate_spectral_entropy`.

        _

    Returns
    -------
    df : pandas.DataFrame
        The dataframe with the spectral entropy column added.
    """
    values = df.apply(lambda x: calculate_spectral_entropy(x[col_peaks], **kwargs), axis=1)
    return df.assign(**{col_result: values})
