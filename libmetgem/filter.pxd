from libcpp.vector cimport vector
from common cimport peak_t

cdef vector[peak_t] filter_data_nogil(float[:,:], float, int, int, int, int) nogil