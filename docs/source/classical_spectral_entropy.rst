================
Spectral entropy
================


Before we calculate spectral entropy, the spectrum need to be centroid first. Which means one fragment ion should only have one peak. When you are focusing on fragment ion's information, the precursor ion may need to be removed from the spectrum before calculating spectral entropy.

The ``calculate_spectral_entropy`` function will do the centroid step and calculate the spectral entropy. For example:

.. code-block:: python

    import numpy as np
    import ms_entropy as me

    peaks = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)

    entropy = me.calculate_spectral_entropy(peaks, clean_spectrum = True, min_ms2_difference_in_da = 0.05)

    print("Spectral entropy is {}.".format(entropy))

If you want to seprate the centroid step and entropy calculation step, you can use ``clean_spectrum`` and ``calculate_spectral_entropy`` function. For example:

.. code-block:: python

    import numpy as np
    import ms_entropy as me

    peaks = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype=np.float32)

    peaks = me.clean_spectrum(peaks, min_ms2_difference_in_da = 0.05)

    entropy = me.calculate_spectral_entropy(peaks, clean_spectrum = False)

    print("Spectral entropy is {}.".format(entropy))


If your spectrum is already centroid, you can skip the ``clean_spectrum`` step. For example:

.. code-block:: python

    import numpy as np
    import scipy.stats

    peaks = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype = np.float32)

    entropy = me.calculate_spectral_entropy(peaks, clean_spectrum = False)

    print("Spectral entropy is {}.".format(entropy))


References
----------

.. autofunction:: ms_entropy.clean_spectrum 
    :noindex:

.. autofunction:: ms_entropy.calculate_spectral_entropy
    :noindex:


