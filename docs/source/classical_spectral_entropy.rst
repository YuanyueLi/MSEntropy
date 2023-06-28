================
Spectral entropy
================


Prior to calculating spectral entropy, the spectrum needs to be centroided, meaning that each fragment ion should only have one corresponding peak. When focusing on fragment ion information, it may be necessary to remove the precursor ion from the spectrum before performing the spectral entropy calculation.

The ``calculate_spectral_entropy`` function carries out the centroiding step and then computes the spectral entropy. Here's an example:

.. code-block:: python

    import numpy as np
    import ms_entropy as me

    peaks = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)

    entropy = me.calculate_spectral_entropy(peaks, clean_spectrum = True, min_ms2_difference_in_da = 0.05)

    print(f"Spectral entropy is {entropy}.")

If you want to separate the centroiding and entropy calculation steps, you can use the ``clean_spectrum`` and ``calculate_spectral_entropy`` functions respectively. Here's how to do it:

.. code-block:: python

    import numpy as np
    import ms_entropy as me

    peaks = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)

    peaks = me.clean_spectrum(peaks, min_ms2_difference_in_da = 0.05)

    entropy = me.calculate_spectral_entropy(peaks, clean_spectrum = False)

    print(f"Spectral entropy is {entropy}.")

If your spectrum is already centroided, you can skip the ``clean_spectrum`` step and directly use the ``entropy`` function from ``scipy.stats`` to compute the spectral entropy. Here's an example:

.. code-block:: python

    import numpy as np
    import scipy.stats

    peaks = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype = np.float32)

    entropy = scipy.stats.entropy(peaks[:, 1])

    print(f"Spectral entropy is {entropy}.")


References
----------

.. autofunction:: ms_entropy.clean_spectrum 
    :noindex:

.. autofunction:: ms_entropy.calculate_spectral_entropy
    :noindex:


