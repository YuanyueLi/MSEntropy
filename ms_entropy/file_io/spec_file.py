#!/usr/bin/env python3
import zipfile
from pathlib import Path
from typing import Union

import numpy as np

from . import lbm2_file, mgf_file, msp_file, mzml_file


def standardize_spectrum(spectrum_dict: dict, standardize_info: dict, keep_all_keys=True):
    """
    Standardize spectrum informat to a standard format provided by standardize_info.

    :param spectrum_dict: spectrum info to standardize, this dict will be modified if keep_all_keys is True.
    :param standardize_info: {wanted_key:[[candidate_keys],default_value,default_type_function]}
    :param keep_all_keys: whether to keep all keys in spectrum_dict. If False, only output wanted keys, else output all keys.

    :return: standardized spectrum info
    """
    if keep_all_keys:
        spec_result = spectrum_dict
    else:
        spec_result = {}
    for key_target, (key_all_candidates, value_default, cast_type) in standardize_info.items():
        key_all_candidates_new = [key_target] + key_all_candidates
        for key_candidate in key_all_candidates_new:
            if key_candidate in spectrum_dict:
                try:
                    if cast_type is None:
                        spec_result[key_target] = spectrum_dict.pop(key_candidate)
                    else:
                        spec_result[key_target] = cast_type(spectrum_dict.pop(key_candidate))
                    break
                except:
                    continue
        else:
            spec_result[key_target] = value_default

        # for key_candidate in key_all_candidates:
        #     spectrum_dict.pop(key_candidate, None)
    return spec_result


def read_one_spectrum(file_input: Union[str, Path],
                      file_type: object = None,
                      **kwargs) -> dict:
    """
    A generator to read one spectrum from file.

    Currently support **.mgf**, **.msp**, **.mzML** and **.lbm2** file.

    The .mgf, .msp can be compressed with .gz, .bz2 or .zip extension.

    The .mgf format is tested with files generated by MSConvert and file downloaded from GNPS.
    The .msp format is tested with files from NIST, files downloaded from MassBank.us and GNPS, and files generated by MS-DIAL.
    The .mzML format is tested with files generated by MSConvert.
    The .lbm2 format is tested with files downloaded from MS-DIAL website.

    :param file_input: The file path of input file, can be str or Path.
    :param file_type: The file type of input file, default is None. Can be "mgf", "msp", "mzml" or "lbm2".
                        If file_type is None, the file type will be determined by file extension.
    :return: a dict contains one spectrum information. The following keys are expected:
                "_ms_level": The MS level of the spectrum, int. Expected value is 1 or 2.
                "_scan_number": The scan number of the spectrum, int. Expected value is 1 or larger.
                "peaks": The peaks of the spectrum, np.array with dtype=np.float32. Expected shape is (N, 2).
    """

    # Determine file type
    if file_type is None:
        file_type = guess_file_type_from_file_name(file_input)
    if file_type is None:
        raise ValueError("Cannot determine file type from file name: {}".format(file_input))

    if file_type == "mgf":
        spectral_generator = mgf_file.read_one_spectrum(file_input, **kwargs)
    elif file_type == "msp":
        spectral_generator = msp_file.read_one_spectrum(file_input, **kwargs)
    elif file_type == "mzml":
        spectral_generator = mzml_file.read_one_spectrum(file_input, **kwargs)
    elif file_type == "lbm2":
        spectral_generator = lbm2_file.read_one_spectrum(file_input, **kwargs)
    else:
        raise ValueError("Unknown file type: {}".format(file_type))

    for spectrum in spectral_generator:
        yield spectrum


def guess_file_type_from_file_name(filename_input):
    file_type = None
    filename_input = str(filename_input)

    # For zip file, only select the first file for the zip file
    if filename_input[-4:].lower() == ".zip":
        fzip_all = zipfile.ZipFile(filename_input)
        fzip_list = zipfile.ZipFile.namelist(fzip_all)
        # Only select the first file for the zip file
        if fzip_list[0][-4:].lower() == ".msp":
            file_type = "msp"
        elif fzip_list[0][-4:].lower() == ".mgf":
            file_type = "mgf"

    # For .gz file
    elif filename_input[-3:].lower() == ".gz":
        if filename_input[-8:-3].lower() == ".mzml":
            file_type = "mzml"
        elif filename_input[-7:-3].lower() == ".msp":
            file_type = "msp"
        elif filename_input[-7:-3].lower() == ".mgf":
            file_type = "mgf"

    # For .bz2 file
    elif filename_input[-4:].lower() == ".bz2":
        if filename_input[-8:-4].lower() == ".msp":
            file_type = "msp"
        elif filename_input[-8:-4].lower() == ".mgf":
            file_type = "mgf"

    else:
        if filename_input[-4:].lower() == ".msp":
            file_type = "msp"
        elif filename_input[-4:].lower() == ".mgf":
            file_type = "mgf"
        elif filename_input[-5:].lower() == ".mzml":
            file_type = "mzml"
        elif filename_input[-5:].lower() == ".hdf5":
            file_type = "hdf5"
        elif filename_input[-4:].lower() == ".raw":
            file_type = "raw"
        elif filename_input[-5:].lower() == ".lbm2":
            file_type = "lbm2"
        elif filename_input.split(".")[-2].lower() == "msp":
            file_type = "msp"
        elif filename_input.split(".")[-2].lower() == "mgf":
            file_type = "mgf"
        elif filename_input.split(".")[-2].lower() == "mzml":
            file_type = "mzml"

    return file_type
