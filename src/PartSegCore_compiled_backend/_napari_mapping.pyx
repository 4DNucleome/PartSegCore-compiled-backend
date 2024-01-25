# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, embedsignature=True

import numpy as np

cimport numpy as cnp

from cython.operator import dereference
from cython.parallel import prange
from libcpp.unordered_map cimport unordered_map

ctypedef fused label_types:
    cnp.uint8_t
    cnp.int8_t
    cnp.uint16_t
    cnp.int16_t
    cnp.uint32_t
    cnp.int32_t
    cnp.uint64_t
    cnp.int64_t

ctypedef fused out_types_mod:
    cnp.uint8_t
    cnp.uint16_t

ctypedef fused out_types:
    out_types_mod
    cnp.float32_t

def _zero_preserving_modulo_seq(label_types[:] labels,  int modulo, int to_zero, out_types_mod[:] out):
    cdef Py_ssize_t i
    cdef Py_ssize_t n = labels.shape[0]
    for i in range(n):
        if labels[i] == to_zero:
            out[i] = 0
        else:
            out[i] = ((modulo + ((labels[i] - 1) % modulo)) % modulo) + 1

def _zero_preserving_modulo_par(label_types[:] labels, int modulo, int to_zero, out_types_mod[:] out):
    cdef Py_ssize_t i
    cdef Py_ssize_t n = labels.shape[0]
    for i in prange(n, nogil=True):
        if labels[i] == to_zero:
            out[i] = 0
        else:
            out[i] = ((modulo + ((labels[i] - 1) % modulo)) % modulo) + 1


def _map_array_par(label_types[:] labels, dict dkt, out_types[:] out, out_types def_val=0):
    # build the map from the input and output vectors
    cdef Py_ssize_t i, n_map, n_array
    cdef unordered_map[label_types, out_types].iterator it
    cdef unordered_map[label_types, out_types] lut
    for key, val in dkt.items():
        if key is not None:
            lut[<label_types>key] = <out_types> val
        else:
            def_val = <out_types> val
    # apply the map to the array
    n_array = labels.shape[0]
    for i in prange(n_array, nogil=True): #
        it = lut.find(labels[i])
        if it != lut.end():
            out[i] = dereference(it).second
        else:
            out[i] = def_val


def _map_array_seq(label_types[:] labels, dict dkt, out_types[:] out, out_types def_val=0):
    # build the map from the input and output vectors
    cdef Py_ssize_t i, n_map, n_array
    cdef unordered_map[label_types, out_types].iterator it
    cdef unordered_map[label_types, out_types] lut
    for key, val in dkt.items():
        if key is not None:
            lut[<label_types>key] = <out_types> val
        else:
            def_val = <out_types> val
    # apply the map to the array
    n_array = labels.shape[0]
    for i in range(n_array): #
        it = lut.find(labels[i])
        if it != lut.end():
            out[i] = dereference(it).second
        else:
            out[i] = def_val

# cdef class MapWrapper:
#     cdef unordered_map[label_types, out_types] lut
