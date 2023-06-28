==================
Entropy similarity
==================


Before initiating entropy similarity calculations, it is critical to first clean the spectrum. Particularly, it is highly recommended to remove any noise ions. Our tests on the NIST20 and Massbank.us databases have shown that eliminating ions with m/z values higher than the precursor ion's m/z - 1.6 considerably enhances the spectral identification performance.

We offer ``calculate_entropy_similarity`` and ``calculate_unweighted_entropy_similarity`` functions to calculate the entropy similarity between two spectra. By default, these functions centroid the spectra and remove noise ions with intensities lower than 1% of the highest intensity ion. If this behavior isn't suited to your needs, parameters can be adjusted; go to the References section for more information.

.. warning:: 
    The ``calculate_entropy_similarity`` and ``calculate_unweighted_entropy_similarity`` functions clean the spectra by default. If you prefer not to clean the spectra, please set ``clean_spectra`` to ``False``. However, when you opt for this, ensure the spectra are pre-centroided, the minimum m/z difference between any two peaks in one spectrum is larger than 2 * ``ms2_tolerance_in_da``, and in every spectrum, the sum of all peaks' intensity is 1. Otherwise, the results may be incorrect.

.. code-block:: python

    import numpy as np
    import ms_entropy as me

    peaks_query = np.array([[69.071, 7.917962], [86.066, 1.021589], [86.0969, 100.0]], dtype = np.float32)
    peaks_reference = np.array([[41.04, 37.16], [69.07, 66.83], [86.1, 999.0]], dtype = np.float32)

    # Calculate unweighted entropy similarity.
    unweighted_similarity = me.calculate_unweighted_entropy_similarity(peaks_query, peaks_reference, ms2_tolerance_in_da = 0.05)
    print(f"Unweighted entropy similarity: {unweighted_similarity}.")

    # Calculate entropy similarity.
    similarity = me.calculate_entropy_similarity(peaks_query, peaks_reference, ms2_tolerance_in_da = 0.05)
    print(f"Entropy similarity: {similarity}.")


References
----------

.. autofunction:: ms_entropy.calculate_entropy_similarity 
    :noindex: 

.. autofunction:: ms_entropy.calculate_unweighted_entropy_similarity  
    :noindex:
    
.. autofunction:: ms_entropy.apply_weight_to_intensity 
    :noindex:
    