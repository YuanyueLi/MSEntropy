#!/usr/bin/env python3
import numpy as np
from .shared import smart_open_file


def read_one_spectrum(file_input, file_encoding='utf-8', **kwargs) -> dict:
    """
    A generator to read one spectrum from .msp file.

    :param filename_input: The file path of input file, can be str or Path.
                            The file can be compressed with .gz, .bz2 or .zip extension.
    :param file_encoding: The encoding of input file, default is 'utf-8'.
    :return: a dict contains one spectrum information.    
    """
    fi = smart_open_file(file_input, file_encoding)

    scan_number = 1
    spectrum_info = {
        "_ms_level": 2,
        "_scan_number": scan_number,
        "peaks": [],
    }
    peak_start = False
    peak_num = 0
    for line in fi:
        if not isinstance(line, str):
            line = line.decode()
        line = line.strip()

        # if scan_number > 10000:
        #     break

        if peak_start:
            items = line.split()
            if len(items) >= 2:
                spectrum_info['peaks'].append([items[0], items[1]])
                peak_num -= 1
            if peak_num == 0:
                # spectrum_info['peaks'] = np.array(spectrum_info['peaks']).astype(np.float32)
                yield spectrum_info
                scan_number += 1
                spectrum_info = {
                    "_ms_level": 2,
                    "_scan_number": scan_number,
                    "peaks": [],
                }
                peak_start = False
        else:
            items = line.split(":", maxsplit=1)
            if len(items) == 2:
                key, value = items
                key = key.strip().lower()
                value = value.strip()
                spectrum_info[key] = value

                if key.lower() == "num peaks":
                    peak_start = True
                    peak_num = int(value)
                    continue

    if len(spectrum_info['peaks']) > 0 or len(spectrum_info) > 3:
        yield spectrum_info
    fi.close()
