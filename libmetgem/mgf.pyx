# distutils: language=c++
# # cython: linetrace=True
# # distutils: define_macros=CYTHON_TRACE_NOGIL=1

cimport cython
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libc.stdlib cimport strtof as std_strtof, strtol
from libc.string cimport strncmp, strncpy, strcpy, strcspn, strlen
from libc.stdio cimport fopen, fclose, fgets, FILE

from .common cimport peak_t, arr_from_vector

DEF MZ = 0
DEF INTENSITY = 1
DEF MAX_LINE_SIZE = 131 # 128 characters + '\r\n' + '\0'

cdef extern from "<string.h>" nogil:
    char *strchr (char *string, int c)

IF WIN32:
    DEF CHARSET = "mbcs"

    cdef extern from "<string.h>" nogil:
        char *strlwr (char *string)
ELSE:
    DEF CHARSET = "UTF-8"

    cdef extern from "<ctype.h>" nogil:
        int tolower(int c)

    cdef char* strlwr(char* string) nogil:
        cdef int i=0

        while string[i] != '\0':
            string[i] = tolower(string[i])
            i += 1

        return string
    
cdef inline float strtof(char* string, char **endptr) nogil:
    cdef char *ptr = NULL
    
    # Allow comma as decimal separator
    ptr = strchr(string, ',')
    if ptr > string:
        string[ptr-string] = '.'
    return std_strtof(string, endptr)
    
cdef void read_data(char line[MAX_LINE_SIZE], vector[peak_t] *peaklist, FILE *fp) nogil:
    """Read peak list from file.
       First line has already been read and has to be passed as first argument.
       peak list is read into `peaklist`, which is cleared as first step."""
       
    cdef:
        peak_t peak
        double value
        char *ptr = NULL
        
    peaklist.clear()
    while True:
        if strncmp(line, 'END IONS', 8) == 0:
            return
        else:
            value = std_strtof(line, &ptr)
            if value > 0:
                peak.mz = value
                peak.intensity = std_strtof(ptr, NULL)
                peaklist.push_back(peak)
                
        if fgets(line, MAX_LINE_SIZE, fp) == NULL:
            return
            
cdef tuple read_entry(FILE * fp, bint ignore_unknown=False):
    """Read a spectrum entry (params and peaklist) from file
    """
    
    cdef:
        dict params = {}
        char *ptr = NULL
        char line[MAX_LINE_SIZE]
        vector[peak_t] peaklist
        int charge
        size_t pos
        char key[32]
        char value[1024]
        
    
    while fgets(line, MAX_LINE_SIZE, fp) != NULL:
        # Ignore blank lines
        if line[0] == '\n' or line[0] == '\r':
            continue
            
        ptr = strchr(line, '=')
        if ptr > line:
            if strncmp(line, 'PEPMASS', 7) == 0:
                params['pepmass'] = strtof(line+8, NULL)
            elif strncmp(line, 'CHARGE', 6) == 0:
                charge = strtol(line+7, &ptr, 10)
                if strncmp(ptr, '-', 1) == 0:
                    charge *= -1
                params['charge'] = charge
            elif strncmp(line, 'RTINSECONDS', 11) == 0:
                params['rtinsecond'] = strtof(line+12, NULL)
            elif strncmp(line, 'MSLEVEL', 7) == 0:
                params['mslevel'] = strtol(line+8, &ptr, 10)
            elif not ignore_unknown:
                pos = ptr - line
                strncpy(key, line, pos)
                key[pos] = '\0'
                key = strlwr(key)
                strcpy(value, line+pos+1)
                if strlen(value) > 0:
                    pos = strcspn(value, '\r\n')
                    value[pos] = '\0'
                    params[key.decode('UTF-8', 'ignore')] = value.decode('UTF-8', 'ignore')
        else:
            if not params or not 'pepmass' in params:
                # If no pepmass found, skip all ions
                while True:
                    if strncmp(line, 'END IONS', 8) == 0:
                        break
                    if fgets(line, MAX_LINE_SIZE, fp) == NULL:
                        break
                return params, np.empty((0,0), dtype=np.float32)
            else:
                read_data(line, &peaklist, fp)
                if peaklist.size() > 0:
                    return params, np.asarray(arr_from_vector(peaklist))
                else:
                    return params, np.empty((0,0), dtype=np.float32)
                return
                
# @cython.binding(True)
def read(str filename, bint ignore_unknown=False):
    cdef:
        bytes fname_bytes
        char *fname
        tuple entry
        char line[MAX_LINE_SIZE]
        FILE *fp

    fname_bytes = filename.encode(CHARSET)
        
    fname = fname_bytes
        
    fp = fopen(fname, 'r')
    if fp == NULL:
        return

    while fgets(line, MAX_LINE_SIZE, fp) != NULL:
        if strncmp(line, 'BEGIN IONS', 10) == 0:
            entry = read_entry(fp, ignore_unknown)
            if entry:
                yield entry
            
    fclose(fp)