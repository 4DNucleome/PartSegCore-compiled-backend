# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, embedsignature=True

from collections.abc import Sequence

import numpy as np

cimport numpy as cnp

from cython.operator cimport preincrement, predecrement, dereference as deref
from libcpp cimport bool
from libcpp.unordered_set cimport unordered_set
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libcpp.set cimport set as cpp_set
from libcpp.algorithm cimport sort


cdef extern from "triangulate.hpp":

    cdef cppclass Event:
        float x
        float y
        int index
        bool is_left
        Event(float x, float y, int index, bool is_left)
        Event()

    cdef cppclass PairHash:
        size_t operator()(pair[int, int] p) const

    cdef cppclass Point:
        float x
        float y
        Point(float x, float y)
        Point()
        operator==(const Point& other) const

    cdef cppclass Segment:
        Point left
        Point right
        Segment(Point left, Point right)
        Segment()

cdef point_eq(Point a, Point b):
    return a.x == b.x and a.y == b.y

cdef bool cmp_event(Event a, Event b):
    if a.y == b.y:
        return a.x < b.x
    return a.y < b.y

cdef bool cmp_event_x(Event a, Event b):
    if a.x == b.x:
        return a.is_left > b.is_left
    return a.x < b.x




cdef inline bool _on_segment(Point p, Point q, Point r):
    if (q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and
        q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y)):
        return True
    return False


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


cdef int _orientation(Point p, Point q, Point r):
    cdef float val
    val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    if val == 0:
        return 0
    return 1 if val > 0 else 2

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


cdef bool _do_intersect(Segment s1, Segment s2):
    cdef Point p1, q1, p2, q2
    cdef int o1, o2, o3, o4
    p1 = s1.left
    q1 = s1.right
    p2 = s2.left
    q2 = s2.right
    if point_eq(p1, p2) or point_eq(p1, q2) or point_eq(q1, p2) or point_eq(q1, q2):
        return False
    o1 = _orientation(p1, q1, p2)
    o2 = _orientation(p1, q1, q2)
    o3 = _orientation(p2, q2, p1)
    o4 = _orientation(p2, q2, q1)

    if o1 != o2 and o3 != o4:
        return True

    if o1 == 0 and _on_segment(p1, p2, q1):
        return True
    if o2 == 0 and _on_segment(p1, q2, q1):
        return True
    if o3 == 0 and _on_segment(p2, p1, q2):
        return True
    if o4 == 0 and _on_segment(p2, q1, q2):
        return True

    return False


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


cdef cpp_set[Event].iterator pred(cpp_set[Event]& s,  cpp_set[Event].iterator it):
    if it == s.begin():
        return s.end()
    return predecrement(it)

cdef cpp_set[Event].iterator succ(cpp_set[Event]& s,  cpp_set[Event].iterator it):
    return preincrement(it)


cdef unordered_set[pair[int, int], PairHash] _find_intersections(vector[Segment] segments):
    cdef unordered_set[pair[int, int], PairHash] intersections
    cdef vector[Event] events
    cdef cpp_set[Event] active, viewed
    cdef Py_ssize_t i, j
    cdef Event current #, next_event, prev_event
    cdef int current_index
    cdef cpp_set[Event].iterator next_, prev, it
    cdef bool flag

    events.reserve(2 * segments.size())
    for i in range(segments.size()):
        events.push_back(Event(segments[i].left.x, segments[i].left.y, i, True))
        events.push_back(Event(segments[i].right.x, segments[i].right.y, i, False))

    sort(events.begin(), events.end(), cmp_event_x)

    for i in range(events.size()):
        current = events[i]
        current_index = current.index
        if current.is_left:
            next_ = active.lower_bound(current)
            prev = pred(active, next_)
            # next_event = dereference(next_)
            # prev_event = dereference(prev)
            flag = False
            if next_ != active.end() and _do_intersect(segments[current_index], segments[deref(next_).index]):
                if current_index < deref(next_).index:
                    intersections.insert(pair[int, int](current_index, deref(next_).index))
                else:
                    intersections.insert(pair[int, int](deref(next_).index, current_index))
            if prev != active.end() and _do_intersect(segments[current_index], segments[deref(prev).index]):
                if current_index < deref(prev).index:
                    intersections.insert(pair[int, int](current_index, deref(prev).index))
                else:
                    intersections.insert(pair[int, int](deref(prev).index, current_index))
            active.insert(current)
            viewed.insert(current)
        else:
            it = active.find(Event(segments[current_index].left.x, segments[current_index].left.y, current_index, True))
            next_ = succ(active, it)
            prev = pred(active, it)
            # next_event = dereference(next_)
            # prev_event = dereference(prev)
            if next_ != active.end() and prev != active.end() and _do_intersect(segments[deref(next_).index], segments[deref(prev).index]):
                if deref(next_).index < deref(prev).index:
                    intersections.insert(pair[int, int](deref(next_).index, deref(prev).index))
                else:
                    intersections.insert(pair[int, int](deref(prev).index, deref(next_).index))
            active.erase(it)
    return intersections


def find_intersections(segments: Sequence[Sequence[Sequence[float]]]) -> list[tuple[int]]:
    """ Find intersections between segments"""
    cdef vector[Segment] segments_vector
    cdef unordered_set[pair[int, int], PairHash] intersections
    cdef pair[int, int] p
    cdef list result = []
    segments_vector.reserve(len(segments))

    for segment in segments:
        if segment[0][0] < segment[1][0]:
            segments_vector.push_back(Segment(Point(segment[0][0], segment[0][1]), Point(segment[1][0], segment[1][1])))
        else:
            segments_vector.push_back(Segment(Point(segment[1][0], segment[1][1]), Point(segment[0][0], segment[0][1])))
    intersections = _find_intersections(segments_vector)
    for p in intersections:
        result.append((p.first, p.second))
    return result
