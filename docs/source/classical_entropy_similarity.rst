==================
Entropy similarity
==================

Entropy similarity
------------------

Before calculate entropy similarity, the spectrum need to be centroid first. Remove the noise ions is highly recommend. Also, base on our test on NIST20 and Massbank.us database, remove ions have m/z higher than precursor ion's m/z - 1.6 will greatly improve the spectral identification performance.

We provide ``calculate_entropy_similarity`` function to calculate two spectral entropy.

.. code-block:: python

    import numpy as np
    import ms_entropy as me

    peaks_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype = np.float32)
    peaks_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype = np.float32)

    # Calculate entropy similarity.
    similarity = me.calculate_entropy_similarity(peaks_query, peaks_reference, ms2_tolerance_in_da = 0.05)
    print("Entropy similarity:{}.".format(similarity))


Unweighted entropy similarity
-----------------------------

Before calculate entropy similarity, the spectrum need to be centroid first. Remove the noise ions is highly recommend. Also, base on our test on NIST20 and Massbank.us database, remove ions have m/z higher than precursor ion's m/z - 1.6 will greatly improve the spectral identification performance.

We provide ``calculate_entropy_similarity`` function to calculate two spectral entropy.

.. code-block:: python

    import numpy as np
    import ms_entropy as me

    peaks_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype = np.float32)
    peaks_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype = np.float32)

    # Calculate entropy similarity.
    similarity = me.calculate_entropy_similarity(peaks_query, peaks_reference, ms2_tolerance_in_da = 0.05)
    print("Entropy similarity:{}.".format(similarity))



References
----------

.. autofunction:: ms_entropy.calculate_entropy_similarity 
    :noindex: 

.. autofunction:: ms_entropy.calculate_unweighted_entropy_similarity  
    :noindex:
    
.. autofunction:: ms_entropy.apply_weight_to_intensity 
    :noindex:
    