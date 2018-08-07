from pyteomics import mgf
from pyteomics.auxiliary import PyteomicsError
from .cosine import MZ, INTENSITY
import numpy as np

def filter_data(data, mz_parent, min_intensity, parent_filter_tolerance, matched_peaks_window,
                min_matched_peaks_search):

    # Filter low mass peaks
    data = data[data[:, MZ] >= 50]

    # Filter peaks close to the parent ion's m/z
    data = data[np.logical_or(data[:, MZ] <= mz_parent - parent_filter_tolerance,
                              data[:, MZ] >= mz_parent + parent_filter_tolerance)]

    if data.size > 0:
        # Keep only peaks higher than threshold
        data = data[data[:, INTENSITY] >= min_intensity * data[:, INTENSITY].max() / 100]

    if data.size > 0:
        # Window rank filter
        data = data[np.argsort(data[:, INTENSITY])]

        if data.size > 0:
            mz_ratios = data[:, MZ]
            mask = np.logical_and(mz_ratios >= mz_ratios[:, None] - matched_peaks_window,
                                  mz_ratios <= mz_ratios[:, None] + matched_peaks_window)
            data = data[np.array([mz_ratios[i] in mz_ratios[mask[i]][-min_matched_peaks_search:]
                                  for i in range(mask.shape[0])])]
            del mask

            if data.size > 0:
                # Use square root of intensities to minimize/maximize effects of high/low intensity peaks
                data[:, INTENSITY] = np.sqrt(data[:, INTENSITY]) * 10

                # Normalize data to norm 1
                data[:, INTENSITY] = data[:, INTENSITY] / np.sqrt(data[:, INTENSITY] @ data[:, INTENSITY])

    return data
