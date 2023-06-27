import numpy as np
import pprint
from ms_entropy import FlashEntropySearch

# This is your library spectra, here the "precursor_mz" and "peaks" are required. The "id" is optional.
spectral_library = [
    {
        "id": "Demo spectrum 1",
        "precursor_mz": 150.0,
        "peaks": [[100.0, 1.0], [101.0, 1.0], [103.0, 1.0]]
    }, {
        "id": "Demo spectrum 2",
        "precursor_mz": 200.0,
        "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32),
        "metadata": "ABC"
    }, {
        "id": "Demo spectrum 3",
        "precursor_mz": 250.0,
        "peaks": np.array([[200.0, 1.0], [101.0, 1.0], [202.0, 1.0]], dtype=np.float32),
        "XXX": "YYY",
    }, {
        "precursor_mz": 350.0,
        "peaks": [[100.0, 1.0], [101.0, 1.0], [302.0, 1.0]]}
]
# This is your query spectrum, here the "peaks" is required, the "precursor_mz" is required for identity search.
query_spectrum = {
    "precursor_mz": 150.0,
    "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32)
}

flash_entropy = FlashEntropySearch()
# Please note that the library spectra will be re-sorted by precursor_mz for fast identity search.
spectral_library = flash_entropy.build_index(spectral_library)

# Clean the spectrum. This step is required!
query_spectrum['peaks'] = flash_entropy.clean_spectrum_for_search(precursor_mz=query_spectrum['precursor_mz'],
                                                                  peaks=query_spectrum['peaks'])

# Perform the identity search.
entropy_similarity = flash_entropy.identity_search(precursor_mz=query_spectrum['precursor_mz'],
                                                   peaks=query_spectrum['peaks'],
                                                   ms1_tolerance_in_da=0.01, ms2_tolerance_in_da=0.02)
# Output the top 3 matches.
topn_match = flash_entropy.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)
print('-' * 20, 'Identity search', '-' * 20)
pprint.pprint(topn_match)

# Perform the open search.
entropy_similarity = flash_entropy.open_search(peaks=query_spectrum['peaks'], ms2_tolerance_in_da=0.02)
# Output the top 3 matches.
topn_match = flash_entropy.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)
print('-' * 20, 'Open search', '-' * 20)
pprint.pprint(topn_match)

# Perform the neutral loss search.
entropy_similarity = flash_entropy.neutral_loss_search(precursor_mz=query_spectrum['precursor_mz'],
                                                       peaks=query_spectrum['peaks'],
                                                       ms2_tolerance_in_da=0.02)
# Output the top 3 matches.
topn_match = flash_entropy.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)
print('-' * 20, 'Neutral loss search', '-' * 20)
pprint.pprint(topn_match)

# Perform the hybrid search.
entropy_similarity = flash_entropy.hybrid_search(precursor_mz=query_spectrum['precursor_mz'], peaks=query_spectrum['peaks'], ms2_tolerance_in_da=0.02)
# Output the top 3 matches.
topn_match = flash_entropy.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)
print('-' * 20, 'Hybrid search', '-' * 20)
pprint.pprint(topn_match)
