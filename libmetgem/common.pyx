# distutils: language=c++

cimport cython
from cython.view cimport array as cvarray
from libcpp.vector cimport vector

DEF MZ = 0
DEF INTENSITY = 1

@cython.boundscheck(False)
@cython.wraparound(False)
cdef float[:,:] arr_from_vector(vector[peak_t] v):
    cdef:
        int i
        float[:,:] a = cvarray(shape=(v.size(), 2), itemsize=sizeof(float), format='f')
        
    for i in range(v.size()):
        a[i, MZ] = v[i].mz
        a[i, INTENSITY] = v[i].intensity
        
    return a