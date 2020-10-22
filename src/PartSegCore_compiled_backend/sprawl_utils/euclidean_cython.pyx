# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, embedsignature=True
# cython: language_level=3


from cpython.mem cimport PyMem_Free
from numpy cimport float64_t, int8_t, uint8_t

from .distance_utils cimport Point, component_types

include "put_borders_in_queue.pyx"


def calculate_euclidean(np.ndarray[uint8_t, ndim=3] object_area, np.ndarray[uint8_t, ndim=3] base_object,
                        np.ndarray[int8_t, ndim=2] neighbourhood, np.ndarray[float64_t, ndim=1] distance,
                        distance_cache=None):
    """
    Calculate euclidean watersheed for one core object

    :param object_area: Area in which euclidean watershhed sholud be calculated
    :param base_object: Core object from which watersheed start
    :param neighbourhood: negihbourhood defined as array of size (neigbourhood_size, 3)
    :param distance: array for distances of negibours. Need have size (neighbourhood_szie).
        Used for handling spacin in image.
    :return: distance from core object
    """
    cdef np.ndarray[uint8_t, ndim=3] consumed_area = np.copy(base_object)
    cdef np.ndarray[float64_t, ndim=3] result, d_cache
    cdef Size x_size, y_size, z_size, array_pos, x, y, z, xx, yy, zz
    cdef Py_ssize_t count = 0
    cdef char neigh_length = neighbourhood.shape[0]
    cdef int neigh_it
    cdef my_queue[Point] current_points
    cdef Point p, p1
    z_size = object_area.shape[0]
    y_size = object_area.shape[1]
    x_size = object_area.shape[2]
    result = np.full((z_size, y_size, x_size), np.inf, dtype=np.float64)
    result[base_object != 0] = 0
    if distance_cache is None:
        d_cache = result
    else:
        d_cache = distance_cache
    put_borders_in_queue(current_points, base_object, neighbourhood)
    while not current_points.empty():
        p = current_points.front()
        current_points.pop()
        count += 1
        """if consumed_area[p.z, p.y, p.x] > 0:
            continue"""
        consumed_area[p.z, p.y, p.x] = 0
        for neigh_it in range(neigh_length):
            z = p.z + neighbourhood[neigh_it, 0]
            y = p.y + neighbourhood[neigh_it, 1]
            x = p.x + neighbourhood[neigh_it, 2]
            if x < 0 or y < 0 or z < 0 or x >= x_size or y >= y_size or z >= z_size:
                continue
            if object_area[z, y, x] == 0:
                continue
            if d_cache[z, y, x] >= result[p.z, p.y, p.x] + distance[neigh_it]:
                if result[z, y, x] > result[p.z, p.y, p.x] + distance[neigh_it]:
                    result[z, y, x] = result[p.z, p.y, p.x] + distance[neigh_it]
                    if consumed_area[z, y, x] == 0:
                        consumed_area[z, y, x] = 1
                        p1.z = z
                        p1.y = y
                        p1.x = x
                        current_points.push(p1)
    return result
