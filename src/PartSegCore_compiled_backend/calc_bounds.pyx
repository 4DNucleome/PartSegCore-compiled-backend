# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, embedsignature=True
"""
This module contains the cython implementation of the ``calc_bounds`` function.

currently supported dimensions: 2, 3, 4, 5
currently supported label types: uint8, uint16, uint32
"""


from typing import Optional
import numpy as np

cimport numpy as cnp

ctypedef fused label_types:
    cnp.uint8_t
    cnp.uint16_t
    cnp.uint32_t

def calc_bounds(labels: np.ndarray, components_num: Optional[int]=None):
    """
    Calculate the bounds of each component.
    wrapper around ``calc_boundsX`` function for different number of dimensions

    Parameters
    ----------
    labels : ndarray
        The labels of the components.
    components_num : int, optional
        The number of components.

    Returns
    -------
    bounds : (ndarray, ndarray)
        The bounds of each component.
    """
    if components_num is None:
        components_num = int(np.max(labels))
    return {
        2: calc_bounds2,
        3: calc_bounds3,
        4: calc_bounds4,
        5: calc_bounds5,
    }[labels.ndim](labels, components_num)


def calc_bounds5(cnp.ndarray[label_types, ndim=5] labels, components_num: Py_ssize_t):
    """
    Calculate the bounds of the specified labels in a 5-dimensional array.

    Parameters
    ----------
    labels: numpy.ndarray[label_types, ndim=5]
        The 5-dimensional array containing the labels.
    components_num: Py_ssize_t
        The number of components to calculate the bounds for.

    Returns
    -------
    bounds: tuple[numpy.ndarray[cnp.int16_t, ndim=2], numpy.ndarray[cnp.int16_t, ndim=2]]
        The bounds of each component.
    """
    cdef Py_ssize_t x, y, z, t, s
    cdef Py_ssize_t x_max = labels.shape[4]
    cdef Py_ssize_t y_max = labels.shape[3]
    cdef Py_ssize_t z_max = labels.shape[2]
    cdef Py_ssize_t t_max = labels.shape[1]
    cdef Py_ssize_t s_max = labels.shape[0]
    cdef label_types label_val
    cdef cnp.ndarray[cnp.int16_t, ndim=2] min_bound = np.full((components_num + 1, 5), max(x_max, y_max, z_max, t_max, s_max) + 5, dtype=np.int16)
    cdef cnp.ndarray[cnp.int16_t, ndim=2] max_bound = np.full((components_num + 1, 5), -1, dtype=np.int16)

    for s in range(0, s_max):
        for t in range(0, t_max):
            for z in range(0, z_max):
                for y in range(0,y_max):
                    for x in range(0,x_max):
                        label_val = labels[s, t, z, y, x]
                        if label_val == 0:
                            continue
                        min_bound[label_val, 0] = min(min_bound[label_val, 0], s)
                        min_bound[label_val, 1] = min(min_bound[label_val, 1], t)
                        min_bound[label_val, 2] = min(min_bound[label_val, 2], z)
                        min_bound[label_val, 3] = min(min_bound[label_val, 3], y)
                        min_bound[label_val, 4] = min(min_bound[label_val, 4], x)

                        max_bound[label_val, 0] = max(max_bound[label_val, 0], s)
                        max_bound[label_val, 1] = max(max_bound[label_val, 1], t)
                        max_bound[label_val, 2] = max(max_bound[label_val, 2], z)
                        max_bound[label_val, 3] = max(max_bound[label_val, 3], y)
                        max_bound[label_val, 4] = max(max_bound[label_val, 4], x)

    return min_bound, max_bound


def calc_bounds4(cnp.ndarray[label_types, ndim=4] labels, components_num: Py_ssize_t):
    """
    Calculate the bounds of the specified labels in a 4-dimensional array.

    Parameters
    ----------
    labels: numpy.ndarray[label_types, ndim=4]
        The 4-dimensional array containing the labels.
    components_num: Py_ssize_t
        The number of components to calculate the bounds for.

    Returns
    -------
    bounds: tuple[numpy.ndarray[cnp.int16_t, ndim=2], numpy.ndarray[cnp.int16_t, ndim=2]]
        The bounds of each component.
    """
    cdef Py_ssize_t x, y, z, t
    cdef Py_ssize_t x_max = labels.shape[3]
    cdef Py_ssize_t y_max = labels.shape[2]
    cdef Py_ssize_t z_max = labels.shape[1]
    cdef Py_ssize_t t_max = labels.shape[0]
    cdef label_types label_val
    cdef cnp.ndarray[cnp.int16_t, ndim=2] min_bound = np.full((components_num + 1, 4), max(x_max, y_max, z_max, t_max) + 5, dtype=np.int16)
    cdef cnp.ndarray[cnp.int16_t, ndim=2] max_bound = np.full((components_num + 1, 4), -1, dtype=np.int16)

    for t in range(0, t_max):
        for z in range(0, z_max):
            for y in range(0,y_max):
                for x in range(0,x_max):
                    label_val = labels[t, z, y, x]
                    if label_val == 0:
                        continue
                    min_bound[label_val, 0] = min(min_bound[label_val, 0], t)
                    min_bound[label_val, 1] = min(min_bound[label_val, 1], z)
                    min_bound[label_val, 2] = min(min_bound[label_val, 2], y)
                    min_bound[label_val, 3] = min(min_bound[label_val, 3], x)

                    max_bound[label_val, 0] = max(max_bound[label_val, 0], t)
                    max_bound[label_val, 1] = max(max_bound[label_val, 1], z)
                    max_bound[label_val, 2] = max(max_bound[label_val, 2], y)
                    max_bound[label_val, 3] = max(max_bound[label_val, 3], x)
    return min_bound, max_bound

def calc_bounds3(cnp.ndarray[label_types, ndim=3] labels, components_num: Py_ssize_t):
    """
    Calculate the bounds of the specified labels in a 3-dimensional array.

    Parameters
    ----------
    labels: numpy.ndarray[label_types, ndim=3]
        The 3-dimensional array containing the labels.
    components_num: Py_ssize_t
        The number of components to calculate the bounds for.

    Returns
    -------
    bounds: tuple[numpy.ndarray[cnp.int16_t, ndim=2], numpy.ndarray[cnp.int16_t, ndim=2]]
        The bounds of each component.
    """
    cdef Py_ssize_t x, y, z
    cdef Py_ssize_t x_max = labels.shape[2]
    cdef Py_ssize_t y_max = labels.shape[1]
    cdef Py_ssize_t z_max = labels.shape[0]
    cdef label_types label_val
    cdef cnp.ndarray[cnp.int16_t, ndim=2] min_bound = np.full((components_num + 1, 3), max(x_max, y_max, z_max) + 5, dtype=np.int16)
    cdef cnp.ndarray[cnp.int16_t, ndim=2] max_bound = np.full((components_num + 1, 3), -1, dtype=np.int16)


    for z in range(0, z_max):
        for y in range(0,y_max):
            for x in range(0,x_max):
                label_val = labels[z, y, x]
                if label_val == 0:
                    continue
                min_bound[label_val, 0] = min(min_bound[label_val, 0], z)
                min_bound[label_val, 1] = min(min_bound[label_val, 1], y)
                min_bound[label_val, 2] = min(min_bound[label_val, 2], x)

                max_bound[label_val, 0] = max(max_bound[label_val, 0], z)
                max_bound[label_val, 1] = max(max_bound[label_val, 1], y)
                max_bound[label_val, 2] = max(max_bound[label_val, 2], x)
    return min_bound, max_bound


def calc_bounds2(cnp.ndarray[label_types, ndim=2] labels, components_num: Py_ssize_t):
    """
    Calculate the bounds of the specified labels in a 2-dimensional array.

    Parameters
    ----------
    labels: numpy.ndarray[label_types, ndim=2]
        The 2-dimensional array containing the labels.
    components_num: Py_ssize_t
        The number of components to calculate the bounds for.

    Returns
    -------
    bounds: tuple[numpy.ndarray[cnp.int16_t, ndim=2], numpy.ndarray[cnp.int16_t, ndim=2]]
        The bounds of each component.
    """
    cdef Py_ssize_t x, y
    cdef Py_ssize_t x_max = labels.shape[1]
    cdef Py_ssize_t y_max = labels.shape[0]
    cdef label_types label_val
    cdef cnp.ndarray[cnp.int16_t, ndim=2] min_bound = np.full((components_num + 1, 2), max(x_max, y_max) + 5, dtype=np.int16)
    cdef cnp.ndarray[cnp.int16_t, ndim=2] max_bound = np.full((components_num + 1, 2), -1, dtype=np.int16)


    for y in range(0,y_max):
        for x in range(0,x_max):
            label_val = labels[y, x]
            if label_val == 0:
                continue
            min_bound[label_val, 0] = min(min_bound[label_val, 0], y)
            min_bound[label_val, 1] = min(min_bound[label_val, 1], x)

            max_bound[label_val, 0] = max(max_bound[label_val, 0], y)
            max_bound[label_val, 1] = max(max_bound[label_val, 1], x)

    return min_bound, max_bound
