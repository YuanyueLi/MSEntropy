import numpy as np
import pprint
from ms_entropy import FlashEntropySearch

# This is your library spectra, here the "precursor_mz" and "peaks" are required. The "id" is optional.
spectral_library = [{
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
    "peaks": [[100.0, 1.0], [101.0, 1.0], [302.0, 1.0]]}]
# This is your query spectrum, here the "peaks" is required, the "precursor_mz" is required for identity search.
query_spectrum = {"precursor_mz": 150.0,
                  "peaks": np.array([[100.0, 1.0], [101.0, 1.0], [102.0, 1.0]], dtype=np.float32)}

entropy_search = FlashEntropySearch()
# Please note that the library spectra will be re-sorted by precursor_mz for fast identity search.
spectral_library = entropy_search.build_index(spectral_library)

# Perform the identity search.
search_result = entropy_search.search(method="identity",
                                      precursor_mz=query_spectrum['precursor_mz'],
                                      peaks=query_spectrum['peaks'])
entropy_similarity = search_result['identity_search']
# Output the top 3 matches.
topn_match = entropy_search.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)
print('-' * 20, 'Identity search', '-' * 20)
pprint.pprint(topn_match)

# Perform the open search.
search_result = entropy_search.search(method="open",
                                      precursor_mz=query_spectrum['precursor_mz'],
                                      peaks=query_spectrum['peaks'])
entropy_similarity = search_result['open_search']
# Output the top 3 matches.
topn_match = entropy_search.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)
print('-' * 20, 'Open search', '-' * 20)
pprint.pprint(topn_match)

# Perform the neutral loss search.
search_result = entropy_search.search(method="neutral_loss",
                                      precursor_mz=query_spectrum['precursor_mz'],
                                      peaks=query_spectrum['peaks'])
entropy_similarity = search_result['neutral_loss_search']
# Output the top 3 matches.
topn_match = entropy_search.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)
print('-' * 20, 'Neutral loss search', '-' * 20)
pprint.pprint(topn_match)

# Perform the hybrid search.
search_result = entropy_search.search(method="hybrid",
                                      precursor_mz=query_spectrum['precursor_mz'],
                                      peaks=query_spectrum['peaks'])
entropy_similarity = search_result['hybrid_search']
# Output the top 3 matches.
topn_match = entropy_search.get_topn_matches(entropy_similarity, topn=3, min_similarity=0.01)
print('-' * 20, 'Hybrid search', '-' * 20)
pprint.pprint(topn_match)
