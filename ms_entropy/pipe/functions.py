import re
import numpy as np
import pandas as pd
from .pipe import pipe_function

__all__ = ["convert_peaks_to_numpy_array", "print_"]


@pipe_function
def convert_peaks_to_numpy_array(df, col_peaks="peaks", delim_intra_peak=":\t", delim_inter_peak=";\n", **kwargs):
    """
    Convert the peaks column from string to numpy array.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe to be converted.

    col_peaks : str
        The column name of the peaks column.

    delim_intra_peak : str
        The delimiter between the mz and intensity of a peak. If multiple characters are provided, all of them will be used as delimiters.

    delim_inter_peak : str
        The delimiter between peaks. If multiple characters are provided, all of them will be used as delimiters.

    Returns
    -------
    df : pandas.DataFrame
        The converted dataframe.
    """
    func_split_intra_peak = lambda x: [y for y in re.split(f"[{delim_intra_peak}]", x)]
    func_split_inter_peak = lambda x: np.array([func_split_intra_peak(y) for y in re.split(f"[{delim_inter_peak}]", x)]).astype(np.float32)
    new_value = df[col_peaks].apply(func_split_inter_peak)
    return df.assign(**{col_peaks: new_value})


@pipe_function
def print_(df, *args, **kwargs):
    print(df, *args, **kwargs)
    return df
