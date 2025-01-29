# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, embedsignature=True

from collections.abc import Sequence

import numpy as np
import cython

cimport numpy as cnp

from libcpp cimport bool
from libcpp.unordered_set cimport unordered_set
from libcpp.utility cimport pair
from libcpp.vector cimport vector


cdef extern from "triangulation/point.hpp" namespace "partsegcore::point":
    cdef cppclass Point:
        float x
        float y
        Point()
        Point(float x, float y)
        bool operator==(const Point& other) const
        bool operator!=(const Point& other) const

    cdef cppclass Vector:
        float x;
        float y;

    cdef cppclass Segment:
        Point bottom
        Point top
        Segment()
        Segment(Point bottom, Point top)


cdef extern from "triangulation/intersection.hpp" namespace "partsegcore::intersection":
    cdef cppclass Event:
        float x
        float y
        int index
        bool is_top
        Event()
        Event(float x, float y, int index, bool is_top)

    cdef cppclass OrderedPair:
        int first
        int second
        OrderedPair()
        OrderedPair(int first, int second)

    bool _on_segment(const Point& p, const Point& q, const Point& r)
    int _orientation(const Point& p, const Point& q, const Point& r)
    bool _do_intersect(const Segment& s1, const Segment& s2)
    unordered_set[OrderedPair] _find_intersections(const vector[Segment]& segments)
    vector[Point] _find_intersection(const Segment& s1, const Segment& s2)


cdef extern from "triangulation/triangulate.hpp" namespace "partsegcore::triangulation":

    cdef cppclass Triangle:
        size_t x
        size_t y
        size_t z
        Triangle()
        Triangle(int x, int y, int z)

    cdef cppclass MonotonePolygon:
        Point top
        Point bottom
        vector[Point] left
        vector[Point] right
        MonotonePolygon()
        MonotonePolygon(Point top, Point bottom, vector[Point] left, vector[Point] right)


    cdef cppclass PointTriangle:
        Point p1
        Point p2
        Point p3

    cdef cppclass PathTriangulation:
        vector[Triangle] triangles
        vector[Point] centers
        vector[Vector] offsets


    bool _is_convex(const vector[Point]& polygon)
    vector[Triangle] _triangle_convex_polygon(const vector[Point]& polygon)
    bool left_to_right(const Segment& s1, const Segment& s2)
    vector[Point] find_intersection_points(const vector[Point]& segments)
    vector[PointTriangle] triangulate_monotone_polygon(const MonotonePolygon& polygon)
    pair[vector[Triangle], vector[Point]] triangulate_polygon_face(const vector[Point]& polygon) except + nogil
    pair[vector[Triangle], vector[Point]] triangulate_polygon_face(const vector[vector[Point]]& polygon_list) except + nogil
    PathTriangulation triangulate_path_edge(const vector[Point]& path, bool closed, float limit, bool bevel) except + nogil
    vector[vector[Point]] split_polygon_on_repeated_edges(const vector[Point]& polygon) except + nogil


ctypedef fused float_types:
    cnp.float32_t
    cnp.float64_t


def on_segment(p: Sequence[float], q: Sequence[float], r: Sequence[float]) -> bool:
    """ Check if point q is on segment pr

    Parameters
    ----------
    p: sequence of 2 floats
        beginning of segment
    q: sequence of 2 floats
        point to check
    r: sequence of 2 floats
        end of segment

    Returns
    -------
    bool:
        True if q is on segment pr
    """
    return _on_segment(
        Point(p[0], p[1]),
        Point(q[0], q[1]),
        Point(r[0], r[1])
        )

def orientation(p: Sequence[float], q: Sequence[float], r: Sequence[float]) -> int:
    """ Check orientation of 3 points

    Parameters
    ----------
    p: sequence of 2 floats
        first point
    q: sequence of 2 floats
        second point
    r: sequence of 2 floats
        third point

    Returns
    -------
    int:
        0 if points are collinear, 1 if clockwise, 2 if counterclockwise
    """
    return _orientation(
        Point(p[0], p[1]),
        Point(q[0], q[1]),
        Point(r[0], r[1])
        )


def do_intersect(s1: Sequence[Sequence[float]], s2: Sequence[Sequence[float]]) -> bool:
    """ Check if two segments intersect

    Parameters
    ----------
    s1: sequence of 2 sequences of 2 floats
        first segment
    s2: sequence of 2 sequences of 2 floats
        second segment

    Returns
    -------
    bool:
        True if segments intersect
    """
    return _do_intersect(
        Segment(Point(s1[0][0], s1[0][1]), Point(s1[1][0], s1[1][1])),
        Segment(Point(s2[0][0], s2[0][1]), Point(s2[1][0], s2[1][1]))
        )


def find_intersections(segments: Sequence[Sequence[Sequence[float]]]) -> list[tuple[int, int]]:
    """ Find intersections between segments"""
    cdef vector[Segment] segments_vector
    cdef unordered_set[OrderedPair] intersections
    cdef OrderedPair p

    segments_vector.reserve(len(segments))
    for segment in segments:
        segments_vector.push_back(Segment(Point(segment[0][0], segment[0][1]), Point(segment[1][0], segment[1][1])))
    intersections = _find_intersections(segments_vector)

    return [(p.first, p.second) for p  in intersections]


def find_intersection_point(s1: Sequence[Sequence[float]], s2: Sequence[Sequence[float]]) -> list[tuple[float, float]]:
    """ Find intersection between two segments

    Parameters
    ----------
    s1: sequence of 2 sequences of 2 floats
        first segment
    s2: sequence of 2 sequences of 2 floats
        second segment

    Returns
    -------
    sequence of 2 floats:
        intersection point
    """
    cdef vector[Point] p_li = _find_intersection(
        Segment(Point(s1[0][0], s1[0][1]), Point(s1[1][0], s1[1][1])),
        Segment(Point(s2[0][0], s2[0][1]), Point(s2[1][0], s2[1][1]))
        )
    return [(p.x, p.y) for p in p_li]


def is_convex(polygon: Sequence[Sequence[float]]) -> bool:
    """ Check if polygon is convex"""
    cdef vector[Point] polygon_vector
    cdef pair[bool, vector[int]] result

    polygon_vector.reserve(len(polygon))
    for point in polygon:
        polygon_vector.push_back(Point(point[0], point[1]))

    return _is_convex(polygon_vector)

def triangle_convex_polygon(polygon: Sequence[Sequence[float]])  -> list[tuple[int, int, int]]:
    cdef vector[Point] polygon_vector
    cdef vector[Triangle] result
    cdef Point p1, p2

    polygon_vector.reserve(len(polygon))
    polygon_vector.push_back(Point(polygon[0][0], polygon[0][1]))
    for point in polygon[1:]:
        p1 = polygon_vector[polygon_vector.size() - 1]
        p2 = Point(point[0], point[1])
        if p1 != p2:
            # prevent from adding polygon edge of width 0
            polygon_vector.push_back(p2)

    result = _triangle_convex_polygon(polygon_vector)
    return [(triangle.x, triangle.y, triangle.z) for triangle in result]



def triangulate_polygon_py(polygon: Sequence[Sequence[float]]) -> tuple[list[tuple[int, int, int]], list[tuple[float, float]]]:
    """ Triangulate polygon"""
    cdef vector[Point] polygon_vector
    cdef Point p1, p2
    cdef pair[vector[Triangle], vector[Point]] result

    polygon_vector.reserve(len(polygon))
    polygon_vector.push_back(Point(polygon[0][0], polygon[0][1]))
    for point in polygon[1:]:
        p1 = polygon_vector[polygon_vector.size() - 1]
        p2 = Point(point[0], point[1])
        if p1 != p2:
            # prevent from adding polygon edge of width 0
            polygon_vector.push_back(p2)

    result = triangulate_polygon_face(polygon_vector)
    return [(triangle.x, triangle.y, triangle.z) for triangle in result.first], [(point.x, point.y) for point in result.second]


def triangulate_polygon_numpy(cnp.ndarray[float_types, ndim=2] polygon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """ Triangulate polygon"""
    cdef vector[Point] polygon_vector
    cdef Point p1, p2
    cdef pair[vector[Triangle], vector[Point]] result

    polygon_vector.reserve(polygon.shape[0])
    polygon_vector.push_back(Point(polygon[0, 0], polygon[0, 1]))
    for point in polygon[1:]:
        p1 = polygon_vector[polygon_vector.size() - 1]
        p2 = Point(point[0], point[1])
        if p1 != p2:
            # prevent from adding polygon edge of width 0
            polygon_vector.push_back(p2)

    result = triangulate_polygon_face(polygon_vector)
    return (
        np.array([(triangle.x, triangle.y, triangle.z) for triangle in result.first], dtype=np.uintp),
        np.array([(point.x, point.y) for point in result.second], dtype=np.float32)
    )


def triangulate_polygon_numpy_li(list polygon_li: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """ Triangulate polygon"""
    cdef vector[Point] polygon_vector
    cdef vector[vector[Point]] polygon_vector_list
    cdef Point p1, p2
    cdef pair[vector[Triangle], vector[Point]] result
    cdef size_t i;
    cdef cnp.ndarray[cnp.float32_t, ndim=2] polygon32
    cdef cnp.ndarray[cnp.float64_t, ndim=2] polygon64

    polygon_vector_list.reserve(len(polygon_li))
    for polygon in polygon_li:
        polygon_vector.clear()
        polygon_vector.reserve(polygon.shape[0])
        polygon_vector.push_back(Point(polygon[0, 0], polygon[0, 1]))

        if polygon.dtype == np.float32:
            polygon32 = polygon
            for i in range(1, polygon.shape[0]):
                p1 = polygon_vector.back()
                p2 = Point(polygon32[i, 0], polygon32[i, 1])
                if p1 != p2:
                    # prevent from adding polygon edge of width 0
                    polygon_vector.push_back(p2)
        else:
            polygon64 = polygon
            for i in range(1, polygon.shape[0]):
                p1 = polygon_vector.back()
                p2 = Point(polygon64[i, 0], polygon64[i, 1])
                if p1 != p2:
                    # prevent from adding polygon edge of width 0
                    polygon_vector.push_back(p2)

        if polygon_vector.size() > 1 and polygon_vector.front() == polygon_vector.back():
            polygon_vector.pop_back()
        polygon_vector_list.push_back(polygon_vector)

    result = triangulate_polygon_face(polygon_vector_list)
    return (
        np.array([(triangle.x, triangle.y, triangle.z) for triangle in result.first], dtype=np.uintp),
        np.array([(point.x, point.y) for point in result.second], dtype=np.float32)
    )

def segment_left_to_right_comparator(s1: Sequence[Sequence[float]], s2: Sequence[Sequence[float]]) -> bool:
    """ Compare segments by bottom point"""
    return left_to_right(
        Segment(Point(s1[0][0], s1[0][1]), Point(s1[1][0], s1[1][1])),
        Segment(Point(s2[0][0], s2[0][1]), Point(s2[1][0], s2[1][1]))
        )


def find_intersection_points_py(polygon: Sequence[Sequence[float]]) -> Sequence[Sequence[float]]:
    """ Find intersection points in polygon"""
    cdef vector[Point] polygon_vector
    cdef vector[Point] result

    polygon_vector.reserve(len(polygon))
    for point in polygon:
        polygon_vector.push_back(Point(point[0], point[1]))

    result = find_intersection_points(polygon_vector)
    return [(point.x, point.y) for point in result]


def triangulate_monotone_polygon_py(top: Sequence[float], bottom: Sequence[float], left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> Sequence[Sequence[Sequence[float]]]:
    """ Triangulate monotone polygon"""
    cdef vector[Point] left_vector, right_vector
    cdef vector[PointTriangle] result
    cdef MonotonePolygon mono_polygon

    left_vector.reserve(len(left))
    for point in left:
        left_vector.push_back(Point(point[0], point[1]))

    right_vector.reserve(len(right))
    for point in right:
        right_vector.push_back(Point(point[0], point[1]))

    mono_polygon = MonotonePolygon(Point(top[0], top[1]), Point(bottom[0], bottom[1]), left_vector, right_vector)
    result = triangulate_monotone_polygon(mono_polygon)
    return [
        [(triangle.p1.x, triangle.p1.y), (triangle.p2.x, triangle.p2.y), (triangle.p3.x, triangle.p3.y)]
        for triangle in result
    ]


def triangulate_path_edge_py(path: Sequence[Sequence[float]], closed: bool=False, limit: float=3.0, bevel: bool=False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Triangulate path"""
    cdef vector[Point] path_vector
    cdef PathTriangulation result
    cdef Point p1, p2
    cdef cnp.ndarray[cnp.uint32_t, ndim=2] triangles

    path_vector.reserve(len(path))
    for point in path:
        path_vector.push_back(Point(point[0], point[1]))
    with cython.nogil:
        result = triangulate_path_edge(path_vector, closed, limit, bevel)

    if result.triangles.size() == 0:
        triangles = np.zeros((0, 3), dtype=np.uint32)
    else:
        triangles = np.array([(triangle.x, triangle.y, triangle.z) for triangle in result.triangles], dtype=np.uint32)

    return (
        np.array([(point.x, point.y) for point in result.centers], dtype=np.float32),
        np.array([(offset.x, offset.y) for offset in result.offsets], dtype=np.float32),
        triangles,
    )

def triangulate_path_edge_numpy(cnp.ndarray[cnp.float32_t, ndim=2] path, bool closed=False, float limit=3.0, bool bevel=False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Triangulate path"""
    cdef vector[Point] path_vector
    cdef PathTriangulation result
    cdef Point p1, p2
    cdef cnp.ndarray[cnp.uint32_t, ndim=2] triangles
    cdef cnp.ndarray[cnp.float32_t, ndim=2] offsets, centers
    cdef size_t i, len_path

    len_path = path.shape[0]

    path_vector.reserve(len_path)
    for i in range(len_path):
        path_vector.push_back(Point(path[i, 0], path[i, 1]))
    with cython.nogil:
        result = triangulate_path_edge(path_vector, closed, limit, bevel)

    triangles = np.empty((result.triangles.size(), 3), dtype=np.uint32)
    centers = np.empty((result.centers.size(), 2), dtype=np.float32)
    offsets = np.empty((result.offsets.size(), 2), dtype=np.float32)

    for i in range(result.triangles.size()):
        triangles[i, 0] = result.triangles[i].x
        triangles[i, 1] = result.triangles[i].y
        triangles[i, 2] = result.triangles[i].z

    for i in range(result.centers.size()):
        centers[i, 0] = result.centers[i].x
        centers[i, 1] = result.centers[i].y

    for i in range(result.offsets.size()):
        offsets[i, 0] = result.offsets[i].x
        offsets[i, 1] = result.offsets[i].y

    return (
        centers,
        offsets,
        triangles,
    )


def triangulate_polygon_with_edge_numpy_li(polygon_li: list[np.ndarray], split_edges: bool=False) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """ Triangulate polygon"""
    cdef vector[Point] polygon_vector
    cdef vector[vector[Point]] polygon_vector_list, edge_split_list
    cdef Point p1, p2
    cdef pair[vector[Triangle], vector[Point]] triangulation_result
    cdef vector[PathTriangulation] edge_result
    cdef cnp.ndarray[cnp.uint32_t, ndim=2] triangles, edge_triangles
    cdef size_t triangle_count = 0

    cdef cnp.ndarray[cnp.float32_t, ndim=2] points, edge_offsets, edges_centers, polygon
    cdef size_t i, j, len_path, edge_triangle_count, edge_center_count, edge_triangle_index, edge_center_index

    polygon_vector_list.reserve(len(polygon_li))
    for i in range(len(polygon_li)):
        polygon = polygon_li[i]

        polygon_vector.reserve(polygon.shape[0])
        polygon_vector.push_back(Point(polygon[0, 0], polygon[0, 1]))

        for j in range(1, polygon.shape[0]):
            p1 = polygon_vector.back()
            p2 = Point(polygon[j, 0], polygon[j, 1])
            if p1 != p2:
                # prevent from adding polygon edge of width 0
                polygon_vector.push_back(p2)
        if polygon_vector.size() > 1 and polygon_vector.front() == polygon_vector.back():
            polygon_vector.pop_back()
        polygon_vector_list.push_back(polygon_vector)
        if split_edges:
            with cython.nogil:
                edge_split_list = split_polygon_on_repeated_edges(polygon_vector)
                for edge_li in edge_split_list:
                    edge_result.push_back(triangulate_path_edge(edge_li, True, 3.0, False))
        else:
            with cython.nogil:
                edge_result.push_back(triangulate_path_edge(polygon_vector, True, 3.0, False))

    with cython.nogil:
        triangulation_result = triangulate_polygon_face(polygon_vector_list)


    triangles = np.empty((triangulation_result.first.size(), 3), dtype=np.uint32)
    for i in range(triangulation_result.first.size()):
        triangles[i, 0] = triangulation_result.first[i].x
        triangles[i, 1] = triangulation_result.first[i].y
        triangles[i, 2] = triangulation_result.first[i].z

    points = np.empty((triangulation_result.second.size(), 2), dtype=np.float32)
    for i in range(triangulation_result.second.size()):
        points[i, 0] = triangulation_result.second[i].x
        points[i, 1] = triangulation_result.second[i].y

    edge_triangle_count = 0
    edge_center_count = 0
    for i in range(edge_result.size()):
        edge_triangle_count += edge_result[i].triangles.size()
        edge_center_count += edge_result[i].centers.size()

    edge_triangles = np.empty((edge_triangle_count, 3), dtype=np.uint32)
    edge_offsets = np.empty((edge_center_count, 2), dtype=np.float32)
    edges_centers = np.empty((edge_center_count, 2), dtype=np.float32)

    edge_triangle_index = 0
    edge_center_index = 0
    for i in range(edge_result.size()):
        for j in range(edge_result[i].triangles.size()):
            edge_triangles[edge_triangle_index, 0] = edge_result[i].triangles[j].x + triangle_count
            edge_triangles[edge_triangle_index, 1] = edge_result[i].triangles[j].y + triangle_count
            edge_triangles[edge_triangle_index, 2] = edge_result[i].triangles[j].z + triangle_count
            edge_triangle_index += 1
        triangle_count += edge_result[i].centers.size()

        for j in range(edge_result[i].centers.size()):
            edges_centers[edge_center_index, 0] = edge_result[i].centers[j].x
            edges_centers[edge_center_index, 1] = edge_result[i].centers[j].y
            edge_offsets[edge_center_index, 0] = edge_result[i].offsets[j].x
            edge_offsets[edge_center_index, 1] = edge_result[i].offsets[j].y
            edge_center_index += 1

    return ((
        triangles,
        points,
    ),
    (
        edges_centers,
        edge_offsets,
        edge_triangles,
    )
    )


def split_polygon_on_repeated_edges_py(polygon: Sequence[Sequence[float]]) -> list[list[tuple[float, float]]]:
    """ Split polygon on repeated edges"""
    cdef vector[Point] polygon_vector
    cdef vector[vector[Point]] result
    cdef Point p1, p2

    polygon_vector.reserve(len(polygon))
    for point in polygon:
        polygon_vector.push_back(Point(point[0], point[1]))

    result = split_polygon_on_repeated_edges(polygon_vector)
    return [
        [(point.x, point.y) for point in polygon_vector]
        for polygon_vector in result
    ]
