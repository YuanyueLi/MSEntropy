#!/usr/bin/env python3
import pickle
import time
import datetime
import urllib.request
from pathlib import Path
import numpy as np
import zipfile

from ms_entropy import FlashEntropySearch


url_mona = r'https://mona.fiehnlab.ucdavis.edu/rest/downloads/retrieve/03d5a22c-c1e1-4101-ac70-9a4eae437ef5'


def main():
    path_data = Path(__file__).parent / 'data'
    file_index = path_data / 'index.pkl'
    all_entropy_search = load_library_index(path_data, file_index)

    ion_mode = 'N'
    entropy_search = all_entropy_search[ion_mode]
    ###########################################################################
    # Randomly select a spectrum from the negative library as the query.
    spec_query = {
        'precursor_mz': 611.2246,
        'peaks': np.array([[205.0351, 78.571429], [206.0434, 61.764706], [380.1302, 100.000000], [423.1333, 34.033613]], dtype=np.float32),
    }
    # Clean the spectrum first.
    spec_query['peaks'] = entropy_search.clean_spectrum_for_search(precursor_mz=spec_query["precursor_mz"], peaks=spec_query['peaks'], noise_threshold=0.01)

    ###########################################################################
    # Identity Search with Flash Entropy Search.
    print('*' * 80)
    print('Identity Search with Flash Entropy Search')
    start = time.time()
    identity_result = entropy_search.identity_search(precursor_mz=spec_query['precursor_mz'],
                                                     peaks=spec_query['peaks'],
                                                     ms1_tolerance_in_da=0.01,
                                                     ms2_tolerance_in_da=0.02)
    print(f'Finished identity search in {time.time() - start:.4f} seconds with {len(identity_result)} results.')
    # Find the top 3 matches.
    print('Top 5 matches:')
    topn = 5
    topn_matches = entropy_search.get_topn_matches(identity_result, topn=topn)
    for i, spec in enumerate(topn_matches):
        print(f'Rank {i+1}: {spec["id"]} with score {spec["entropy_similarity"]:.4f}')

    ###########################################################################
    # Open Search with Flash Entropy Search.
    print('*' * 80)
    print('Open Search with Flash Entropy Search')
    start = time.time()
    open_result = entropy_search.open_search(peaks=spec_query['peaks'], ms2_tolerance_in_da=0.02)
    print(f'Finished open search in {time.time() - start:.4f} seconds with {len(open_result)} results.')
    # Find the top 5 matches.
    print('Top 5 matches:')
    topn_matches = entropy_search.get_topn_matches(open_result, topn=topn)
    for i, spec in enumerate(topn_matches):
        print(f'Rank {i+1}: {spec["id"]} with score {spec["entropy_similarity"]:.4f}')

    ###########################################################################
    # Neutral Loss Search with Flash Entropy Search.
    print('*' * 80)
    print('Neutral Loss Search with Flash Entropy Search')
    start = time.time()
    neutral_loss_result = entropy_search.neutral_loss_search(precursor_mz=spec_query['precursor_mz'],
                                                             peaks=spec_query['peaks'],
                                                             ms2_tolerance_in_da=0.02)
    print(f'Finished neutral loss search in {time.time() - start:.4f} seconds with {len(neutral_loss_result)} results.')
    # Find the top 5 matches.
    print('Top 5 matches:')
    topn_matches = entropy_search.get_topn_matches(neutral_loss_result, topn=topn)
    for i, spec in enumerate(topn_matches):
        print(f'Rank {i+1}: {spec["id"]} with score {spec["entropy_similarity"]:.4f}')

    ###########################################################################
    # Hybrid Search with Flash Entropy Search.
    print('*' * 80)
    print('Hybrid Search with Flash Entropy Search')
    start = time.time()
    hybrid_result = entropy_search.hybrid_search(precursor_mz=spec_query['precursor_mz'],
                                                 peaks=spec_query['peaks'],
                                                 ms2_tolerance_in_da=0.02)
    print(f'Finished hybrid search in {time.time() - start:.4f} seconds with {len(hybrid_result)} results.')
    # Find the top 5 matches.
    print('Top 5 matches:')
    topn_matches = entropy_search.get_topn_matches(hybrid_result, topn=topn)
    for i, spec in enumerate(topn_matches):
        print(f'Rank {i+1}: {spec["id"]} with score {spec["entropy_similarity"]:.4f}')


def load_library_index(path_data, file_index):
    # If the index does not exist, build it.
    if not file_index.exists():
        # Download MassBank.us spectra.
        file_mona = download_mona(path_data)

        # Load spectra from NoNA data
        print(f'Loading spectra from {file_mona}, this may take a while.')
        all_spectra = reading_spectra_from_mona(str(file_mona))

        # Build the index.
        print('Building index, this will only need to be done once.')
        all_entropy_search = {}
        for ion_mode in all_spectra:
            print(f'Building index for spectra with ion mode {ion_mode}')
            spectra = all_spectra[ion_mode]
            entropy_search = FlashEntropySearch()  # Create a new FlashEntropySearch object.
            spectra = entropy_search.build_index(spectra)  # Build the index.
            entropy_search.library_id = [spec['id'] for spec in spectra]  # Record the library IDs.
            all_entropy_search[ion_mode] = entropy_search

        # Save the index.
        print('Saving index')
        with open(file_index, 'wb') as f:
            pickle.dump(all_entropy_search, f)
    else:
        # Load the index.
        print('Loading index')
        with open(file_index, 'rb') as f:
            all_entropy_search = pickle.load(f)
    return all_entropy_search


def download_mona(path_data):
    # Download MassBank.us spectra.
    path_data.mkdir(parents=True, exist_ok=True)
    file_mona = path_data / f'mona-{datetime.date.today()}.zip'

    if not file_mona.exists():
        print(f'Downloading {url_mona} to {file_mona}')
        with urllib.request.urlopen(url_mona) as response:
            data = response.read()
            with open(file_mona, 'wb') as f:
                f.write(data)
    return file_mona


def reading_spectra_from_mona(file_mona):
    all_spectra = {}
    for i, spec in enumerate(read_one_spectrum(file_mona)):
        if i % 1000 == 0:
            print(f'Total read {i} spectra', end='\r')
        ion_mode = spec.get("ion_mode", "")
        if ion_mode in {'P', 'N'}:
            try:
                spec['precursormz'] = float(spec['precursormz'])
            except:
                continue

            if ion_mode not in all_spectra:
                all_spectra[ion_mode] = []
            # Record the spectrum.
            all_spectra[ion_mode].append({'peaks': spec['peaks'], 'precursor_mz': spec['precursormz'], 'id': spec['db#']})

    print(f'Loaded {len(all_spectra["P"])} positive spectra and {len(all_spectra["N"])} negative spectra.')
    return all_spectra


def read_one_spectrum(filename_input: str) -> dict:
    """
    Read one spectrum from .msp file.
    :param filename_input: a stream for input.
    :return: a dict contains one spectrum information.
    """
    # Deal with compressed file
    if filename_input.lower()[-4:] == ".zip":
        fzip_all = zipfile.ZipFile(filename_input)
        fzip_list = zipfile.ZipFile.namelist(fzip_all)
        fi = fzip_all.open(fzip_list[0], "r")
    else:
        fi = open(filename_input, "rt", encoding='utf-8')

    _scan_number = 1
    spectrum_info = {
        "_ms_level": 2,
        "_scan_number": _scan_number,
        'peaks': []
    }
    is_adding_information = True
    is_spec_end = False
    for line in fi:
        if not isinstance(line, str):
            line = line.decode()

        line = line.strip()

        if is_adding_information:
            item, _ = _parse_information(line, spectrum_info)
            if item.startswith("num") and item.lower() == "num peaks":
                spectrum_info[item] = int(spectrum_info[item])
                is_adding_information = False
                peak_num = spectrum_info[item]
                if peak_num == 0:
                    is_spec_end = True
        else:
            peaks = spectrum_info['peaks']
            items = line.split()
            peaks.append([items[0], items[1]])
            if len(peaks) == peak_num:
                is_spec_end = True

        if is_spec_end:
            spectrum_info['peaks'] = np.array(peaks).astype(np.float32)
            yield spectrum_info

            # Preparing adding the next one.
            is_adding_information = True
            is_spec_end = False
            _scan_number += 1
            spectrum_info = {
                "_scan_number": _scan_number,
                "_ms_level": 2,
                'peaks': []
            }

    fi.close()


def _parse_information(line: bytes, spec: dict):
    """
    Parse the line in .msp file, update information in the spec
    :param line: The input line.
    :param spec: The output information collection.
    :return: The entry of the added information.
    """
    if line:
        lines = line.split(":")
        if len(lines) > 2:
            item = lines[0]
            cont = ":".join(lines[1:])
        elif len(lines) == 2:
            item, cont = lines
        else:
            return "", ""

        item = item.strip().lower()
        cont = cont.strip()
        spec[item] = cont
        return item, cont
    else:
        return "", ""


if __name__ == '__main__':
    main()
