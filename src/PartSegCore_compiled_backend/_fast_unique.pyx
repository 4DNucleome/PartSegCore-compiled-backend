# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, embedsignature=True

import numpy as np

cimport numpy as np
cimport openmp

from cython.parallel import prange, parallel, threadid

from libcpp.unordered_set cimport unordered_set
from libcpp.vector cimport vector


ctypedef fused numpy_int_types:
    np.uint8_t
    np.uint16_t
    np.uint32_t
    np.int8_t
    np.int16_t
    np.int32_t
    np.int64_t


def unique2d(np.ndarray[numpy_int_types, ndim=2] array):
    cdef unordered_set[numpy_int_types] unique_values
    cdef unordered_set[numpy_int_types] unique_values_local
    cdef vector[unordered_set[numpy_int_types]] unique_values_vector
    cdef Py_ssize_t i, j, n_rows, n_cols, id_num
    cdef numpy_int_types prev_val, v
    n_rows = array.shape[0]
    n_cols = array.shape[1]

    with nogil, parallel():
        unique_values_local = unordered_set[numpy_int_types]()
        id_num = threadid()
        if id_num == 0:
            unique_values_vector.resize(openmp.omp_get_num_threads())

        for i in prange(n_rows):
            prev_val = array[i, 0]
            unique_values_local.insert(prev_val)
            for j in range(n_cols):
                if array[i, j] != prev_val:
                    prev_val = array[i, j]
                    unique_values_local.insert(prev_val)

        unique_values_vector[id_num] = unique_values_local

    for s in unique_values_vector:
        for v in s:
            unique_values.insert(v)
    return np.sort(np.fromiter(unique_values, dtype=array.dtype))


def unique3d(np.ndarray[numpy_int_types, ndim=3] array):
    cdef unordered_set[numpy_int_types] unique_values
    cdef unordered_set[numpy_int_types] unique_values_local
    cdef vector[unordered_set[numpy_int_types]] unique_values_vector
    cdef Py_ssize_t i, j, k, n_rows, n_cols, n_depth, id_num
    cdef numpy_int_types prev_val, v
    n_rows = array.shape[0]
    n_cols = array.shape[1]
    n_depth = array.shape[2]

    with nogil, parallel():
        unique_values_local = unordered_set[numpy_int_types]()
        id_num = threadid()
        if id_num == 0:
            unique_values_vector.resize(openmp.omp_get_num_threads())

        for i in prange(n_rows):
            prev_val = array[i, 0, 0]
            unique_values_local.insert(prev_val)
            for j in range(n_cols):
                for k in range(n_depth):
                    if array[i, j, k] != prev_val:
                        prev_val = array[i, j, k]
                        unique_values_local.insert(prev_val)

        unique_values_vector[id_num] = unique_values_local

    for s in unique_values_vector:
        for v in s:
            unique_values.insert(v)

    return np.sort(np.fromiter(unique_values, dtype=array.dtype))
