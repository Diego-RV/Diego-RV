from collections import defaultdict
import numpy as np

# Function to find the best average MS2 value for MS1 peaks with multiple linked MS2 values - V3
def average_msms(linked_msp, top_n, mz_merge_thresh):
    # Initialize a dictionary to hold the MS2 spectra for each MS1 peak
    peak_msms_data = defaultdict(list)  # group MS2 spectra by their corresponding MS1 peak

    for annotation in linked_msp:
        msms_id = annotation['msms_id']
        peak_id = annotation['peak_id']
        centroided_ms2_peaks = np.array(annotation['centroided_ms2_peaks'], dtype=np.float32)
        raw_ms2_peaks = np.array(annotation['raw_ms2_peaks'], dtype=np.float32)
        msp_peaks = annotation['msp_peaks']

        peak_msms_data[peak_id].append({
            'cluster_id': annotation['cluster_id'],
            'file_id': annotation['file_id'],
            'feature_id': annotation['feature_id'],
            'feature_peak_ids': annotation['feature_peak_ids'],
            'peak_id': peak_id,
            'msms_id': msms_id,
            'average_ms1_mz': annotation['average_ms1_mz'],
            'precursor_mz': annotation['precursor_mz'],
            'precursor_intensity': annotation['precursor_intensity'],
            'precursor_rt': annotation['precursor_rt'],
            'precursor_volume': annotation['precursor_volume'],
            'normalized_area': annotation['normalized_area'],
            'charge_state': annotation['charge_state'],
            'centroided_ms2_peaks': centroided_ms2_peaks,
            'unweighted_entropy_similarity': annotation['unweighted_entropy_similarity'],
            'entropy_similarity': annotation['entropy_similarity'],
            'dot_product': annotation['dot_product'],
            'matches': annotation['matches'],
            'raw_ms2_peaks': raw_ms2_peaks,
            'msp_peaks': msp_peaks
        })

    # Initialize a list to hold the results
    averaged_msms = []

    # Process each MS1 peak that has multiple associated MS2 spectra
    for peak_id, msms_list in peak_msms_data.items():
        if len(msms_list) == 1:
            msms = msms_list[0]
            cluster_id = msms['cluster_id']
            feature_id = msms['feature_id']
            msp_peaks_arrays = msms['msp_peaks']
            entropy_similarity = msms['entropy_similarity']
            avg_centroided_ms2_peaks = msms['centroided_ms2_peaks']
            avg_dot_product = msms['dot_product']
            matches = msms['matches']
            msp_peaks = msms['msp_peaks']
            raw_ms2_peaks = msms['raw_ms2_peaks']

            averaged_msms.append({
                'cluster_id': cluster_id,
                'feature_id': feature_id,
                'peak_id': peak_id,
                'avg_dot_product': avg_dot_product,
                'num_msms': 1,
                'avg_entropy_similarity': entropy_similarity,
                'matches': matches,
                'avg_centroided_ms2_peaks': avg_centroided_ms2_peaks,
                'msp_peaks': msp_peaks,
                'avg_ms2_raw_peaks': raw_ms2_peaks
            })

        else:
            unsorted_msms = [m for m in msms_list if m['entropy_similarity'] is not None]
            if not unsorted_msms:
                continue

            # Sort and take top-N by similarity
            sorted_msms = sorted(unsorted_msms, key=lambda x: x['entropy_similarity'], reverse=True)
            top_msms = sorted_msms[:top_n]

            all_raw_ms2_peaks = np.concatenate([np.array(m['raw_ms2_peaks']) for m in top_msms], axis=0)
            sorted_all_raw_peaks = all_raw_ms2_peaks[all_raw_ms2_peaks[:, 0].argsort()]

            all_cent_peaks = np.concatenate([np.array(m['centroided_ms2_peaks']) for m in top_msms], axis=0)
            sorted_all_cent_peaks = all_cent_peaks[all_cent_peaks[:, 0].argsort()]

            groups_current = []
            current = [sorted_all_cent_peaks[0]]
            for mz_i, intensity_i in sorted_all_cent_peaks[1:]:
                if abs(mz_i - current[-1][0]) <= mz_merge_thresh:
                    current.append([mz_i, intensity_i])
                else:
                    groups_current.append(np.array(current))
                    current = [[mz_i, intensity_i]]
                    groups_current.append(np.array(current))

            # Compute average values and intensities for each group
            averaged_results = []
            for group in groups_current:
            # group[:, 0] → all mz values, group[:, 1] → all intensities
                avg_mz = np.mean(group[:, 0])
                avg_intensity = np.mean(group[:, 1])
                # Append as (mz, intensity) tuple
                averaged_results.append((avg_mz, avg_intensity))
                avg_ms2_cent_peaks = np.array(averaged_results, dtype=float)

            msp_peaks_arrays = [m['msp_peaks'] for m in top_msms]
            cluster_id = [m['cluster_id'] for m in top_msms]
            feature_id = [m['feature_id'] for m in top_msms]
            matches = [m['matches'] for m in top_msms]

            # Convert to the desired format: np.array([[m/z, intensity], ...])
            peaks_query = np.array(sorted_all_raw_peaks, dtype=np.float32)
            peaks_reference = np.array(msp_peaks_arrays[0], dtype=np.float32)
            if peaks_reference.ndim == 3 and peaks_reference.shape[0] == 1:
                peaks_reference = peaks_reference[0]

            # Calculate similarity scores (MS_entropy)
            avg_entropy_similarity = me.calculate_entropy_similarity(peaks_query, peaks_reference)

            # Extract m/z and intensity values
            reference_mz, reference_intensity = zip(*msp_peaks_arrays[0])
            ms2_cent_intensity, ms2_mz_values = zip(*avg_ms2_cent_peaks)

            reference_mz = np.array(reference_mz)
            reference_intensity = np.array(reference_intensity, dtype=np.float32)
            sample_mz_values = np.sort(np.array(ms2_mz_values))

            matched_indices = np.searchsorted(reference_mz, sample_mz_values)
            matched_indices = np.clip(matched_indices, 0, len(reference_mz) - 1)

            sample_spectrum = np.array(ms2_cent_intensity[:len(matched_indices)])
            reference_spectrum = reference_intensity[matched_indices]

            sample_spectrum_norm = sample_spectrum / np.linalg.norm(sample_spectrum)
            reference_spectrum_norm = reference_spectrum / np.linalg.norm(reference_spectrum)

            avg_dot_product = np.dot(sample_spectrum_norm, reference_spectrum_norm)

            averaged_msms.append({
                'cluster_id': cluster_id,
                'feature_id': feature_id,
                'peak_id': peak_id,
                'total_num_msms': len(msms_list),
                'avg_dot_product': avg_dot_product,
                'avg_entropy_similarity': avg_entropy_similarity,
                'avg_centroided_ms2_peaks': avg_ms2_cent_peaks,
                'avg_ms2_raw_peaks' : sorted_all_raw_peaks,
                'msp_peaks': msp_peaks_arrays[0],
                'matches': matches
            })

    return averaged_msms
