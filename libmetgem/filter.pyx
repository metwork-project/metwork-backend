# distutils: language=c++

cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.algorithm cimport sort
from .common cimport peak_t, arr_from_vector

DEF MZ = 0
DEF INTENSITY = 1

cdef bool compareByIntensity(const peak_t &a, const peak_t &b) nogil:
    return a.intensity > b.intensity
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef vector[peak_t] filter_data_nogil(float[:,:] data, float mz_parent, int min_intensity, int parent_filter_tolerance, int matched_peaks_window, int min_matched_peaks_search) nogil:
    cdef Py_ssize_t size = data.shape[0]
    cdef int i, j, count=0
    cdef float mz, intensity, max_intensity=0
    cdef float abs_min_intensity
    cdef vector[peak_t] peaks
    cdef vector[peak_t] peaks2
    cdef peak_t peak
    
    # Sort data array by decreasing intensities
    sort(<peak_t*>&data[0,0], (<peak_t*>&data[0,0]) + size, compareByIntensity)

    # Filter out peaks with mz below 50 Da or with mz in `mz_parent` +- `parent_filter_tolerance` or with intensity < `min_intensity` % of maximum intensity
    # Maximum intensity is calculated from peaks not filtered out by mz filters
    peaks.reserve(size)
    for i in range(size):
        mz = data[i, MZ]
        if 50 <= mz <= mz_parent - parent_filter_tolerance or mz >= mz_parent + parent_filter_tolerance:  # mz filter
            intensity = data[i, INTENSITY]
            if intensity < abs_min_intensity:  # intensity filter
                break
            elif intensity > max_intensity:
                max_intensity = intensity
                abs_min_intensity = min_intensity * max_intensity / 100
            
            peak.mz = mz
            peak.intensity = intensity
            peaks.push_back(peak)
    
    # Window rank filter: For each peak, keep it only if it is in the top `min_matched_peaks_search` peaks in the +/- `matched_peaks_window` range
    size = peaks.size()
    peaks2.reserve(size)
    for i in range(size):
        mz = peaks[i].mz
        count = 0
        for j in range(size):
            if mz - matched_peaks_window <= peaks[j].mz <= mz + matched_peaks_window:
                if j == i:
                    peaks[i].intensity = sqrt(peaks[i].intensity) * 10  # Use square root of intensities to minimize/maximize effects of high/low intensity peaks
                    peaks2.push_back(peaks[i])
                    break
                count += 1
                if count >= min_matched_peaks_search:
                    break
        
    return peaks2 #<float[:peaks.size(),:2]>(<float*>peaks.data())
        
def filter_data(np.ndarray[np.float32_t, ndim=2] data, float mz_parent, int min_intensity, int parent_filter_tolerance, int matched_peaks_window, int min_matched_peaks_search):
    cdef np.ndarray[np.float32_t, ndim=2] filtered
    cdef vector[peak_t] peaks = filter_data_nogil(data, mz_parent, min_intensity, parent_filter_tolerance, matched_peaks_window, min_matched_peaks_search)
    
    if peaks.size() == 0:
        return np.empty((0,2), dtype=np.float32)
        
    filtered = np.asarray(arr_from_vector(peaks))
    
    # Normalize data to norm 1
    filtered[:, INTENSITY] = filtered[:, INTENSITY] / np.sqrt(filtered[:, INTENSITY] @ filtered[:, INTENSITY])
    
    return filtered