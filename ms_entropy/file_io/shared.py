#!/usr/bin/env python3
import gzip
import bz2
import zipfile
import codecs


def smart_open_file(file_input, file_encoding='utf-8'):
    """
    Open file with different compression format.

    :param file_input: The file path of input file, can be str or Path.
                        The file can be compressed with .gz, .bz2 or .zip extension.    
    :param file_encoding: The encoding of input file, default is 'utf-8'.
    :return: The file object.
    """

    file_input = str(file_input)
    if file_input.lower()[-3:] == ".gz":
        fi = gzip.open(file_input)
    elif file_input.lower()[-4:] == ".bz2":
        fi = bz2.open(file_input, "rt")
    elif file_input.lower()[-4:] == ".zip":
        fzip_all = zipfile.ZipFile(file_input)
        fzip_list = zipfile.ZipFile.namelist(fzip_all)
        fi = fzip_all.open(fzip_list[0], "r")
    else:
        fi = codecs.open(file_input, "r", encoding=file_encoding, errors='ignore')
    return fi
