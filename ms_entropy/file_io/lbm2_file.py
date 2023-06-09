#!/usr/bin/env python3
import numpy as np


def read_one_spectrum(file_input, **kwargs) -> dict:
    """
    Read all spectra from .lbm2 file.

    :param file_input: The file path of .lbm2 file.
    :return: A list of dict contains all spectra information.
    """
    try:
        import lz4.block
        import msgpack
    except ImportError:
        raise ImportError("Please install lz4 and msgpack package via `pip install lz4 msgpack`.")

    MSGPAK_INT32_LENGTH = 5
    with open(file_input, "rb") as fi:
        buffer = fi.read()
        buffer_data = msgpack.unpackb(buffer)
        lbm2_data = buffer_data.data

        lz4_data_length = msgpack.unpackb(lbm2_data[:MSGPAK_INT32_LENGTH])
        raw_buffer = lbm2_data[MSGPAK_INT32_LENGTH:]
        decompressed_buffer = lz4.block.decompress(raw_buffer, uncompressed_size=lz4_data_length)

        # print("Unpack all spectra")
        all_raw_spectra = msgpack.unpackb(decompressed_buffer)

    # Construct all spectra
    for i, spec_raw in enumerate(all_raw_spectra):
        adduct = spec_raw[13]  # AdductIonBean
        charge = {0: 1, 1: -1}[spec_raw[6]]  # IonMode
        charge = charge*adduct[3]  # ChargeNumber

        peaks = spec_raw[11]  # MzIntensityCommentBean
        peaks = (np.array(peaks)[:, 1:3]).astype(np.float32)

        spec = {
            "_ms_level": 2,
            "_scan_number": i+1, # ScanNumber starts from 1
            "id": spec_raw[0],  # Id
            "compound_class": spec_raw[2],  # CompoundClass
            "name": spec_raw[3],  # Name
            "rt": spec_raw[4],  # RetentionTime
            "smiles": spec_raw[7],  # Smiles
            "inchikey": spec_raw[8],  # InchiKey
            "precursor_mz": spec_raw[9],  # PrecursorMz
            "peaks": peaks,
            "charge": charge,
            "adduct": adduct[2],
        }
        yield spec
    return all_raw_spectra


if __name__ == "__main__":
    file_input = "/p/FastEntropySearch/next_version/io/test.lbm2"
    all_spectra = [x for x in read_one_spectrum(file_input)]
