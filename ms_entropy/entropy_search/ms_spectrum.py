import numpy as np
import struct

# SPECTRUM_STRUCT_STR = r"QffhH"
SPECTRUM_STRUCT_STR = "<QffhH"
HEADER_SIZE = struct.calcsize(SPECTRUM_STRUCT_STR)  # use instead of 20


def read_spectrum_from_file_stream(stream, offset_end=None):
    """
    Read spectra sequentially from a binary file stream.

    This function is a generator that reads MS/MS spectra one by one from the current position of a binary stream. 
    Each spectrum consists of a fixed-size header followed by a variable-length array of peaks.

    Parameters
    ----------
    stream : file-like object

    offset_end : int or None, optional
        Byte offset indicating the end position in the stream up to which spectra should be read. 
        If ``None`` (default), reading continues until end-of-file is reached.

    Yields
    ------
    dict
        A dictionary describing a single spectrum.

    """
    while offset_end is None or stream.tell() < offset_end:
        data = stream.read(HEADER_SIZE)
        if len(data) < HEADER_SIZE:
            break

        scan, precursor_mz, rt, charge, peak_len = struct.unpack(SPECTRUM_STRUCT_STR, data)
        if precursor_mz < 0:
            raise ValueError(f"Bad precursor_mz at {stream.tell()}, Debug info: {precursor_mz}, {scan}, {rt}, {charge}, {peak_len}.")
        try:
            peaks = np.frombuffer(stream.read(peak_len * 8), dtype=np.float32).reshape((peak_len, 2))
        except Exception:
            raise ValueError(f"Error reading peaks for scan {scan} at position {stream.tell()}")

        yield {
            "scan": scan,
            "precursor_mz": precursor_mz,
            "charge": charge,
            "rt": rt,
            "peaks": peaks,
        }


def write_spectrum_to_file_stream(spectrum, stream):

    """
    Write a single spectrum to a binary file stream.

    This function serializes a spectrum dictionary into its binary representation and writes it to the current position of a binary stream.

    Parameters
    ----------
    spectrum : dict
        Spectrum data to be written. The dictionary is expected to follow the internal spectrum schema used by the library and must be compatible with ``convert_spectrum_to_bytes``.

    stream : file-like object

    Returns
    -------
    int
        Status code indicating success. Always returns ``0``.
    """
    stream.write(convert_spectrum_to_bytes(spectrum))
    return 0


def convert_spectrum_to_bytes(spectrum: dict):
    """
    Convert a spectrum dictionary into its binary representation.

    This function serializes a spectrum into a byte sequence consisting of a fixed-size header followed by a contiguous array of peak data. 
    The resulting bytes object can be written directly to a binary file stream.

    Parameters
    ----------
    spectrum : dict
        Spectrum data to be serialized. 

    Returns
    -------
    bytes
        Binary representation of the spectrum.

    """
    n = spectrum["peaks"].shape[0]
    if n > 0xFFFF:
        raise ValueError("Too many peaks for uint16 length field")

    info = struct.pack(
        SPECTRUM_STRUCT_STR,
        spectrum["scan"],
        spectrum["precursor_mz"],
        spectrum.get("rt", -1),
        spectrum["charge"],
        spectrum["peaks"].shape[0],
    )
    return info + spectrum["peaks"].astype(np.float32).tobytes(order="C")


def convert_bytes_to_spectrum(data: bytes):
    """
    Deserialize a binary spectrum representation into a spectrum dictionary.

    This function converts a byte sequence produced by ``convert_spectrum_to_bytes`` back into an in-memory spectrum representation. 
    The input must contain exactly one serialized spectrum.

    Parameters
    ----------
    data : bytes
        Binary data containing a single serialized spectrum. 
        The byte layout must consist of a fixed-size header followed by a contiguous block of peak data.

    Returns
    -------
    dict
        Spectrum dictionary.

    """
    if len(data) < HEADER_SIZE:
        raise ValueError("Data too short to contain a valid spectrum")

    scan, precursor_mz, rt, charge, peak_len = struct.unpack(SPECTRUM_STRUCT_STR, data[:HEADER_SIZE])
    expected_size = HEADER_SIZE + peak_len * 8
    if len(data) != expected_size:
        raise ValueError(f"Data size {len(data)} does not match expected size {expected_size} for {peak_len} peaks")

    peaks = np.frombuffer(data[HEADER_SIZE:], dtype=np.float32).reshape((peak_len, 2))

    return {
        "scan": scan,
        "precursor_mz": precursor_mz,
        "charge": charge,
        "rt": rt,
        "peaks": peaks,
    }
