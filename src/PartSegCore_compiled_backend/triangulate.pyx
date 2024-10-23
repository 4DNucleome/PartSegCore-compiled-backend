# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, embedsignature=True

from collections.abc import Sequence

import numpy as np

cimport numpy as cnp

from libcpp cimport bool
from libcpp.unordered_set cimport unordered_set
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.algorithm cimport sort


cdef extern from "triangulation/point.hpp" namespace "partsegcore::point":
    cdef cppclass Point:
        float x
        float y
        Point()
        Point(float x, float y)
        bool operator==(const Point& other) const
        bool operator!=(const Point& other) const

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
    Point _find_intersection(const Segment& s1, const Segment& s2)


cdef extern from "triangulation/triangulate.hpp" namespace "partsegcore::triangulation":

    cdef cppclass Triangle:
        int x
        int y
        int z
        Triangle()
        Triangle(int x, int y, int z)

    bool _is_convex(const vector[Point]& polygon)
    vector[Triangle] _triangle_convex_polygon(const vector[Point]& polygon)
    bool left_to_right(const Segment& s1, const Segment& s2)



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


def find_intersection_point(s1: Sequence[Sequence[float]], s2: Sequence[Sequence[float]]) -> tuple[float, float]:
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
    cdef Point p = _find_intersection(
        Segment(Point(s1[0][0], s1[0][1]), Point(s1[1][0], s1[1][1])),
        Segment(Point(s2[0][0], s2[0][1]), Point(s2[1][0], s2[1][1]))
        )
    return (p.x, p.y)


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

cdef vector[Triangle] _triangulate_polygon(vector[Point] polygon):
    cdef vector[Segment] edges, edges_with_intersections
    cdef Py_ssize_t i, j, edges_count
    cdef unordered_set[OrderedPair] intersections
    cdef unordered_map[int, vector[Point]] intersections_points
    cdef pair[int, vector[Point]] p_it
    cdef vector[Triangle] triangles
    cdef vector[Point] intersections_points_vector
    cdef Point p_int
    cdef OrderedPair p

    if _is_convex(polygon):
        return _triangle_convex_polygon(polygon)

    edges.reserve(polygon.size())
    for i in range(polygon.size() - 1):
        edges.push_back(Segment(polygon[i], polygon[i+1]))

    intersections = _find_intersections(edges)
    intersections_points.reserve(intersections.size())
    for p in intersections:
        p_int = _find_intersection(edges[p.first], edges[p.second])
        intersections_points[p.first].push_back(p_int)
        intersections_points[p.second].push_back(p_int)

    edges_count = edges.size()
    for p_it in intersections_points:
        edges_count += p_it.second.size() - 1

    edges_with_intersections.reserve(edges_count)

    for i in range(edges.size()):
        if not intersections_points.count(i):
            edges_with_intersections.push_back(edges[i])
        else:
            intersections_points_vector = intersections_points.at(i)
            intersections_points_vector.push_back(edges[i].bottom)
            intersections_points_vector.push_back(edges[i].top)
            sort(intersections_points_vector.begin(), intersections_points_vector.end())
            for j in range(intersections_points_vector.size() - 1):
                edges_with_intersections.push_back(Segment(intersections_points_vector[j], intersections_points_vector[j+1]))

    return triangles


def triangulate_polygon(polygon: Sequence[Sequence[float]]) -> list[tuple[int, int, int]]:
    """ Triangulate polygon"""
    cdef vector[Point] polygon_vector
    cdef Point p1, p2
    cdef vector[Triangle] result


    polygon_vector.reserve(len(polygon))
    polygon_vector.push_back(Point(polygon[0][0], polygon[0][1]))
    for point in polygon[1:]:
        p1 = polygon_vector[polygon_vector.size() - 1]
        p2 = Point(point[0], point[1])
        if p1 != p2:
            # prevent from adding polygon edge of width 0
            polygon_vector.push_back(p2)

    result = _triangulate_polygon(polygon_vector)
    return [(triangle.x, triangle.y, triangle.z) for triangle in result]


def segment_left_to_right_comparator(s1: Sequence[Sequence[float]], s2: Sequence[Sequence[float]]) -> bool:
    """ Compare segments by bottom point"""
    return left_to_right(
        Segment(Point(s1[0][0], s1[0][1]), Point(s1[1][0], s1[1][1])),
        Segment(Point(s2[0][0], s2[0][1]), Point(s2[1][0], s2[1][1]))
        )
