use std::cmp::Ordering;

type FloatSpec = f32;

fn print_spectrum(info: &str, spectrum_2d: &Vec<[f32; 2]>) {
    println!("{}", info);
    for i in 0..spectrum_2d.len() {
        println!("{}\t{}\t{}", i, spectrum_2d[i][0], spectrum_2d[i][1]);
    }
}

fn compare_by_mz_with_zero_intensity(a: &[FloatSpec; 2], b: &[FloatSpec; 2]) -> Ordering {
    if a[1] > 0.0 && b[1] <= 0.0 {
        return Ordering::Less;
    } else if a[1] <= 0.0 && b[1] > 0.0 {
        return Ordering::Greater;
    }

    if a[0] < b[0] {
        return Ordering::Less;
    } else if a[0] > b[0] {
        return Ordering::Greater;
    }

    return Ordering::Equal;
}

fn sort_spectrum_by_mz_and_zero_intensity(spectrum_2d: &mut Vec<[f32; 2]>) {
    spectrum_2d.sort_by(compare_by_mz_with_zero_intensity);

    let mut spectrum_len = spectrum_2d.len();
    while spectrum_len > 0 && spectrum_2d[spectrum_len - 1][1] <= 0.0 {
        spectrum_len -= 1;
    }
    spectrum_2d.truncate(spectrum_len);
}

fn calculate_spectrum_argsort(spectrum_2d: &Vec<[f32; 2]>) -> Vec<usize> {
    let mut spectrum_argsort: Vec<usize> = (0..spectrum_2d.len()).collect();
    spectrum_argsort
        .sort_unstable_by(|&i, &j| spectrum_2d[j][1].partial_cmp(&spectrum_2d[i][1]).unwrap());
    spectrum_argsort
}

// Check if centroiding is needed
fn need_centroid(
    spectrum_2d: &Vec<[FloatSpec; 2]>,
    min_ms2_difference_in_da: FloatSpec,
    min_ms2_difference_in_ppm: FloatSpec,
) -> bool {
    let mut min_ms2_difference_in_da = min_ms2_difference_in_da;
    for i in 0..(spectrum_2d.len() - 1) {
        if min_ms2_difference_in_ppm > 0.0 {
            min_ms2_difference_in_da = spectrum_2d[i + 1][0] * min_ms2_difference_in_ppm * 1e-6;
        }
        if spectrum_2d[i + 1][0] - spectrum_2d[i][0] < min_ms2_difference_in_da {
            return true;
        }
    }
    false
}
// Centroid the spectrum
fn centroid_spectrum(
    mut spectrum_2d: Vec<[FloatSpec; 2]>,
    min_ms2_difference_in_da: FloatSpec,
    min_ms2_difference_in_ppm: FloatSpec,
) -> Vec<[FloatSpec; 2]> {
    // Calculate the argsort of the spectrum by intensity.
    let spectrum_argsort: Vec<usize> = calculate_spectrum_argsort(&spectrum_2d);

    let mut mz_delta_allowed_left = min_ms2_difference_in_da;
    let mut mz_delta_allowed_right = min_ms2_difference_in_da;

    for &idx in &spectrum_argsort {
        if min_ms2_difference_in_ppm > 0.0 {
            mz_delta_allowed_left = spectrum_2d[idx][0] * min_ms2_difference_in_ppm * 1e-6;
            mz_delta_allowed_right = spectrum_2d[idx][0] / (1.0 - min_ms2_difference_in_ppm * 1e-6);
        }
        if spectrum_2d[idx][1] > 0.0 {
            // Find left board for current peak
            let mut idx_left = if idx > 0 { idx - 1 } else { idx };
            while idx_left > 0
                && spectrum_2d[idx][0] - spectrum_2d[idx_left][0] <= mz_delta_allowed_left
            {
                idx_left -= 1;
            }

            // Find right board for current peak
            let mut idx_right = idx + 1;
            while idx_right < spectrum_2d.len()
                && spectrum_2d[idx_right][0] - spectrum_2d[idx][0] <= mz_delta_allowed_right
            {
                idx_right += 1;
            }

            // Merge the peaks in the board
            let mut intensity_sum = 0.0;
            let mut intensity_weighted_sum = 0.0;
            for i in (idx_left + 1)..idx_right {
                intensity_sum += spectrum_2d[i][1];
                intensity_weighted_sum += spectrum_2d[i][1] * spectrum_2d[i][0];
                spectrum_2d[i][1] = 0.0;
            }

            // Write the new peak into the output spectrum
            spectrum_2d[idx][0] = intensity_weighted_sum / intensity_sum;
            spectrum_2d[idx][1] = intensity_sum;
        }
    }

    // Sort by MZ and filter out zero intensity
    sort_spectrum_by_mz_and_zero_intensity(&mut spectrum_2d);

    spectrum_2d
}

fn clean_spectrum(
    spectrum: &mut [FloatSpec],
    min_mz: FloatSpec,
    max_mz: FloatSpec,
    noise_threshold: FloatSpec,
    min_ms2_difference_in_da: FloatSpec,
    min_ms2_difference_in_ppm: FloatSpec,
    max_peak_num: usize,
    normalize_intensity: bool,
) -> Vec<[FloatSpec; 2]> {
    let mut spectrum: Vec<[FloatSpec; 2]> = spectrum
        .chunks_exact(2)
        .map(|chunk| [chunk[0], chunk[1]])
        .collect();
    let spectrum_length = spectrum.len();
    let min_mz = if min_mz < 0.0 { 0.0 } else { min_mz };

    for data in spectrum.iter_mut() {
        if data[0] <= min_mz || (max_mz > 0.0 && data[0] >= max_mz) {
            data[1] = 0.0;
        }
    }
    sort_spectrum_by_mz_and_zero_intensity(&mut spectrum);

    while need_centroid(
        &spectrum,
        min_ms2_difference_in_da,
        min_ms2_difference_in_ppm,
    ) {
        spectrum = centroid_spectrum(
            spectrum,
            min_ms2_difference_in_da,
            min_ms2_difference_in_ppm,
        );
    }

    if noise_threshold > 0.0 {
        let max_intensity = spectrum
            .iter()
            .map(|&data| data[1])
            .fold(FloatSpec::NEG_INFINITY, FloatSpec::max);
        let noise_threshold_intensity = noise_threshold * max_intensity;
        for data in spectrum.iter_mut() {
            if data[1] < noise_threshold_intensity {
                data[1] = 0.0;
            }
        }
    }

    if max_peak_num > 0 && max_peak_num < spectrum_length {
        // Sort the spectrum by intensity and keep only top `max_peak_num` elements
        spectrum.sort_unstable_by(|a, b| b[1].partial_cmp(&a[1]).unwrap());
        spectrum.truncate(max_peak_num);
    }
    sort_spectrum_by_mz_and_zero_intensity(&mut spectrum);

    if normalize_intensity {
        let sum_intensity: FloatSpec = spectrum.iter().map(|&data| data[1]).sum();
        if sum_intensity > 0.0 {
            for data in spectrum.iter_mut() {
                data[1] /= sum_intensity;
            }
        }
    }

    spectrum
}

fn main() {
    let mut spectrum_2d: Vec<[f32; 2]> = vec![
        [41.04, 0.3716],
        [0., 0.3716],
        [69.070, 7.917962],
        [69.070, -7.917962],
        [69.071, 100.],
        [86.0969, 66.83],
        [86.01, 10.],
    ];
    println!("Origin spectrum:");
    print_spectrum("Origin:", &mut spectrum_2d);

    let mut spectrum_data = spectrum_2d.concat();
    let spectrum_2d_final = clean_spectrum(&mut spectrum_data, 0.0, 0.0, 0.01, 0.02, 0.0, 0, true);

    println!("Final spectrum:");
    print_spectrum("Final:", &spectrum_2d_final);
}
