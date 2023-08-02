test_that("msentropy_similarity works", {
    x <- cbind(mz = c(169.071, 186.066, 186.0769),
               intensity = c(7.917962, 1.021589, 100.0))
    y <- cbind(mz = c(120.212, 169.071, 186.066),
               intensity <- c(37.16, 66.83, 999.0))
    res <- msentropy_similarity(x, y)
    ref <- calculate_entropy_similarity(x, y, ms2_tolerance_in_da = 0.02,
                                        ms2_tolerance_in_ppm = -1,
                                        clean_spectra = TRUE,
                                        min_mz = 0, max_mz = 1000,
                                        noise_threshold = 0.01,
                                        max_peak_num = 100)
    expect_equal(res, ref)
    res <- msentropy_similarity(x, y, ms2_tolerance_in_da = -1,
                                ms2_tolerance_in_ppm = 20)
    ref <- calculate_entropy_similarity(x, y, ms2_tolerance_in_da = -1,
                                        ms2_tolerance_in_ppm = 20,
                                        clean_spectra = TRUE,
                                        min_mz = 0, max_mz = 1000,
                                        noise_threshold = 0.01,
                                        max_peak_num = 100)
    expect_equal(res, ref)

    res <- msentropy_similarity(x, y, weighted = FALSE)
    ref <- calculate_unweighted_entropy_similarity(
        x, y, ms2_tolerance_in_da = 0.02,
        ms2_tolerance_in_ppm = -1,
        clean_spectra = TRUE,
        min_mz = 0, max_mz = 1000,
        noise_threshold = 0.01,
        max_peak_num = 100)
    expect_equal(res, ref)
})
