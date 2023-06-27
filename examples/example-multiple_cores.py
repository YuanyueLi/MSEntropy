import numpy as np
from ms_entropy import FlashEntropySearch
import multiprocessing as mp

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

query_spectra_list = [query_spectrum] * 100

entropy_search = FlashEntropySearch()
# Please note that the library spectra will be re-sorted by precursor_mz for fast identity search.
spectral_library = entropy_search.build_index(spectral_library)
entropy_search.save_memory_for_multiprocessing()

THREADS = 4


def func_max_similarity(precursor_mz, peaks):
    entropy_search = func_max_similarity.entropy_search
    similarity = entropy_search.search(precursor_mz=precursor_mz, peaks=peaks, method="neutral_loss")
    return np.max(similarity["neutral_loss_search"])


def init_worker(entropy_search, ):
    func_max_similarity.entropy_search = entropy_search
    return None


pool = mp.Pool(THREADS, initializer=init_worker, initargs=(entropy_search, ))
max_entropy_similarity = pool.starmap(func_max_similarity, [(spectrum["precursor_mz"], spectrum["peaks"]) for spectrum in query_spectra_list])

print(max_entropy_similarity)
