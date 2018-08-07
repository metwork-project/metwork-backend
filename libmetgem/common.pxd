from libcpp.vector cimport vector

ctypedef struct peak_t:
    float mz
    float intensity
    
cdef float[:,:] arr_from_vector(vector[peak_t])