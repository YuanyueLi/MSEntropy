use std::cmp::Ordering;
use wasm_bindgen::prelude::*;

type FloatSpec = f32;

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

fn sort_spectrum_by_mz_and_zero_intensity(spectrum_2d: &mut Vec<[FloatSpec; 2]>) {
    spectrum_2d.sort_by(compare_by_mz_with_zero_intensity);

    let mut spectrum_len = spectrum_2d.len();
    while spectrum_len > 0 && spectrum_2d[spectrum_len - 1][1] <= 0.0 {
        spectrum_len -= 1;
    }
    spectrum_2d.truncate(spectrum_len);
}

fn calculate_spectrum_argsort(spectrum_2d: &Vec<[f32; 2]>) -> Vec<usize> {
    let mut spectrum_argsort: Vec<usize> = (0..spectrum_2d.len()).collect();
    spectrum_argsort.sort_unstable_by(|&i, &j| spectrum_2d[j][1].partial_cmp(&spectrum_2d[i][1]).unwrap());
    spectrum_argsort
}

// Check if centroiding is needed
fn need_centroid(spectrum_2d: &Vec<[FloatSpec; 2]>, min_ms2_difference_in_da: FloatSpec, min_ms2_difference_in_ppm: FloatSpec) -> bool {
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
fn centroid_spectrum(spectrum_2d: &mut Vec<[FloatSpec; 2]>, min_ms2_difference_in_da: FloatSpec, min_ms2_difference_in_ppm: FloatSpec) {
    // Calculate the argsort of the spectrum by intensity.
    let spectrum_argsort: Vec<usize> = calculate_spectrum_argsort(&spectrum_2d);

    let mut mz_delta_allowed_left = min_ms2_difference_in_da;
    let mut mz_delta_allowed_right = min_ms2_difference_in_da;

    for idx in spectrum_argsort {
        if min_ms2_difference_in_ppm > 0.0 {
            mz_delta_allowed_left = spectrum_2d[idx][0] * min_ms2_difference_in_ppm * 1e-6;
            mz_delta_allowed_right = spectrum_2d[idx][0] / (1.0 - min_ms2_difference_in_ppm * 1e-6);
        }
        if spectrum_2d[idx][1] > 0.0 {
            // Find left board for current peak
            let mut idx_left = if idx > 0 { idx - 1 } else { idx };
            while idx_left > 0 && spectrum_2d[idx][0] - spectrum_2d[idx_left][0] <= mz_delta_allowed_left {
                idx_left -= 1;
            }

            // Find right board for current peak
            let mut idx_right = idx + 1;
            while idx_right < spectrum_2d.len() && spectrum_2d[idx_right][0] - spectrum_2d[idx][0] <= mz_delta_allowed_right {
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
    sort_spectrum_by_mz_and_zero_intensity(spectrum_2d);
}

fn clean_spectrum(
    spectrum: &mut Vec<[FloatSpec; 2]>,
    min_mz: FloatSpec,
    max_mz: FloatSpec,
    noise_threshold: FloatSpec,
    min_ms2_difference_in_da: FloatSpec,
    min_ms2_difference_in_ppm: FloatSpec,
    max_peak_num: usize,
    normalize_intensity: bool,
) {
    let spectrum_length = spectrum.len();
    let min_mz = if min_mz < 0.0 { 0.0 } else { min_mz };

    for data in spectrum.iter_mut() {
        if data[0] <= min_mz || (max_mz > 0.0 && data[0] >= max_mz) {
            data[1] = 0.0;
        }
    }
    sort_spectrum_by_mz_and_zero_intensity(spectrum);

    while need_centroid(&spectrum, min_ms2_difference_in_da, min_ms2_difference_in_ppm) {
        centroid_spectrum(spectrum, min_ms2_difference_in_da, min_ms2_difference_in_ppm);
    }

    if noise_threshold > 0.0 {
        let max_intensity = spectrum.iter().map(|&data| data[1]).fold(0., FloatSpec::max);
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
    sort_spectrum_by_mz_and_zero_intensity(spectrum);

    if normalize_intensity {
        let sum_intensity: FloatSpec = spectrum.iter().map(|&data| data[1]).sum();
        if sum_intensity > 0.0 {
            for data in spectrum.iter_mut() {
                data[1] /= sum_intensity;
            }
        }
    }
}

fn calculate_unweighted_entropy_similarity(
    peaks_a: &mut Vec<[FloatSpec; 2]>,
    peaks_b: &mut Vec<[FloatSpec; 2]>,
    ms2_tolerance_in_da: f32,
    ms2_tolerance_in_ppm: f32,
    clean_spectra: bool,
    min_mz: f32,
    max_mz: f32,
    noise_threshold: f32,
    max_peak_num: usize,
) -> f32 {
    if clean_spectra {
        clean_spectrum(
            peaks_a,
            min_mz,
            max_mz,
            noise_threshold,
            2.0 * ms2_tolerance_in_da,
            2.0 * ms2_tolerance_in_ppm,
            max_peak_num,
            true,
        );
        clean_spectrum(
            peaks_b,
            min_mz,
            max_mz,
            noise_threshold,
            2.0 * ms2_tolerance_in_da,
            2.0 * ms2_tolerance_in_ppm,
            max_peak_num,
            true,
        );
    }

    if peaks_a.is_empty() || peaks_b.is_empty() {
        return 0.0;
    }

    let mut a = 0;
    let mut b = 0;
    let mut similarity = 0f32;

    while a < peaks_a.len() && b < peaks_b.len() {
        let mass_delta_da = peaks_a[a][0] - peaks_b[b][0];
        let ms2_tolerance_in_da = if ms2_tolerance_in_ppm > 0.0 {
            ms2_tolerance_in_ppm * peaks_a[a][0] * 1e-6
        } else {
            ms2_tolerance_in_da
        };

        if mass_delta_da < -ms2_tolerance_in_da {
            a += 1;
        } else if mass_delta_da > ms2_tolerance_in_da {
            b += 1;
        } else {
            let peak_a_intensity = peaks_a[a][1];
            let peak_b_intensity = peaks_b[b][1];
            let peak_ab_intensity = peak_a_intensity + peak_b_intensity;
            similarity +=
                peak_ab_intensity * peak_ab_intensity.log2() - peak_a_intensity * peak_a_intensity.log2() - peak_b_intensity * peak_b_intensity.log2();
            a += 1;
            b += 1;
        }
    }
    similarity = similarity / 2.0;
    if similarity < 0.0 {
        similarity = 0.0;
    } else if similarity > 1.0 {
        similarity = 1.0;
    }
    return similarity;
}

#[wasm_bindgen]
extern {
    pub fn alert(s: &str);
}

#[wasm_bindgen]
pub fn greet(name: &str) {
    alert(&format!("Hello, {}!", name));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let spectrum = [
            [41.04, 0.3716],
            [0., 0.3716],
            [69.070, 7.917962],
            [69.070, -7.917962],
            [69.071, 100.],
            [86.0969, 66.83],
            [86.01, 10.],
        ];
        let mut spectrum_2d: Vec<[FloatSpec; 2]> = spectrum.concat().chunks_exact(2).map(|chunk| [chunk[0], chunk[1]]).collect();

        clean_spectrum(&mut spectrum_2d, 0.0, 0.0, 0.01, 0.02, 0.0, 0, true);
        println!("{:?}", spectrum_2d);
    }

    #[test]
    fn entropy_similarity() {
        let spectrum = [
            [41.04, 0.3716],
            [0., 0.3716],
            [69.070, 7.917962],
            [69.070, -7.917962],
            [69.071, 100.],
            [86.0969, 66.83],
            [86.01, 10.],
        ];
        let mut spectrum_2d_a: Vec<[FloatSpec; 2]> = spectrum.concat().chunks_exact(2).map(|chunk| [chunk[0], chunk[1]]).collect();
        let mut spectrum_2d_b: Vec<[FloatSpec; 2]> = spectrum_2d_a.iter().map(|&data| [data[0], data[1]]).collect();
        let similarity = calculate_unweighted_entropy_similarity(&mut spectrum_2d_a, &mut spectrum_2d_b, 0.02, -1f32, true, -1f32, -1f32, 0.01, 0);

        println!("{}", similarity);
    }
}
