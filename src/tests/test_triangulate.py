import pytest

from PartSegCore_compiled_backend.triangulate import (
    on_segment,
    orientation,
    do_intersect,
    find_intersections,
    find_intersection_point,
    is_convex,
    triangulate_polygon,
    triangle_convex_polygon,
)


def test_on_segment():
    assert on_segment((0, 0), (0, 1), (0, 2))
    assert not on_segment((0, 0), (1, 1), (0, 3))


def test_orientation():
    assert orientation((0, 0), (0, 1), (0, 2)) == 0
    assert orientation((0, 0), (0, 2), (0, 1)) == 0
    assert orientation((0, 2), (0, 0), (0, 1)) == 0
    assert orientation((0, 0), (0, 1), (1, 2)) == 1
    assert orientation((0, 0), (0, 1), (-1, 2)) == 2


def test_do_intersect():
    assert do_intersect(((0, -1), (0, 1)), ((1, 0), (-1, 0)))
    assert not do_intersect(((0, -1), (0, 1)), ((1, 0), (2, 0)))


def test_find_intersections():
    r"""
    First test case:
    (1, 0) --- (1, 1)
    |           |
    (0, 0) --- (0, 1)

    Second test case:
    (1, 0) --- (1, 1)
        \     /
         \   /
          \ /
           X
          / \
         /   \
        /     \
    (0, 0) --- (0, 1)
    """
    assert find_intersections([[(0, 0), (0, 1)], [(0, 1), (1, 1)], [(1, 1), (1, 0)], [(1, 0), (0, 0)]]) == []
    assert find_intersections([[(0, 0), (0, 1)], [(0, 1), (1, 0)], [(1, 0), (1, 1)], [(1, 1), (0, 0), (0, 1)]]) == [
        (1, 3)
    ]


@pytest.mark.parametrize(
    'segments, expected',
    [
        # No intersections, simple square
        ([[(0, 0), (0, 1)], [(0, 1), (1, 1)], [(1, 1), (1, 0)], [(1, 0), (0, 0)]], []),
        # One intersection, crossing diagonals
        ([[(0, 0), (2, 2)], [(2, 0), (0, 2)]], [(0, 1)]),
        # Multiple intersections, complex shape
        (
            [[(0, 0), (2, 2)], [(2, 0), (0, 2)], [(1, 0), (1, 2)], [(0, 1), (2, 1)]],
            {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)},
        ),
        # No intersections, non-intersecting lines
        ([[(0, 0), (1, 1)], [(2, 2), (3, 3)]], []),
        # One intersection, T-shaped intersection
        ([[(0, 0), (2, 0)], [(1, -1), (1, 1)]], [(0, 1)]),
        # Multiple intersections, grid shape
        (
            [
                [(0, 0), (2, 0)],
                [(0, 1), (2, 1)],
                [(0, 2), (2, 2)],
                [(0, 0), (0, 2)],
                [(1, 0), (1, 2)],
                [(2, 0), (2, 2)],
            ],
            {(0, 4), (1, 3), (1, 4), (1, 5), (2, 4)},
        ),
    ],
    ids=[
        'No intersections, simple square',
        'One intersection, crossing diagonals',
        'Multiple intersections, complex shape',
        'No intersections, non-intersecting lines',
        'One intersection, T-shaped intersection',
        'Multiple intersections, grid shape',
    ],
)
def test_find_intersections_param(segments, expected):
    assert set(find_intersections(segments)) == set(expected)


# def test_find_intersection_point():
#     assert find_intersection_point(((0, 0), (2, 2)), ((0, 2), (2, 0))) == (1, 1)


@pytest.mark.parametrize(
    'segment1, segment2, expected',
    [
        (((0, 0), (2, 2)), ((0, 2), (2, 0)), (1, 1)),  # Intersecting diagonals
        (((0, 0), (1, 1)), ((1, 0), (0, 1)), (0.5, 0.5)),  # Intersecting diagonals
        (((0, 0), (2, 2)), ((2, 2), (0, 2)), (2, 2)),  # Touching at one point
        (((0, 0), (2, 2)), ((1, 0), (1, 2)), (1, 1)),  # Vertical line intersection
        (((0, 1), (2, 1)), ((2, 0), (2, 2)), (2, 1)),  # Touching in the middle
    ],
    ids=[
        'Intersecting diagonals (1, 1)',
        'Intersecting diagonals (0.5, 0.5)',
        'Touching at one point (2, 2)',
        'Vertical line intersection (1, 1)',
        'Touching in the middle (2, 1)',
    ],
)
def test_find_intersection_point(segment1, segment2, expected):
    assert find_intersection_point(segment1, segment2) == expected


def test_is_convex():
    assert is_convex([(0, 0), (0, 1), (1, 1), (1, 0)])
    assert not is_convex([(0, 0), (0, 1), (1, 0), (1, 1)])

    assert is_convex([(0, 0), (0, 0.5), (0, 1), (1, 1), (1, 0)])
    assert not is_convex([(0, 0), (0.1, 0.5), (0, 1), (1, 1), (1, 0)])


def test_triangle_convex_polygon():
    assert triangle_convex_polygon([(0, 0), (0, 1), (1, 1), (1, 0)]) == [(0, 1, 2), (0, 2, 3)]
    assert triangle_convex_polygon([(0, 0), (0, 0.5), (0, 1), (1, 1), (1, 0)]) == [(0, 2, 3), (0, 3, 4)]


def test_triangulate_polygon():
    assert triangulate_polygon([(0, 0), (0, 1), (1, 1), (1, 0)]) == [(0, 1, 2), (0, 2, 3)]
    assert triangulate_polygon([(0, 0), (0, 10), (1, 10), (1, 0)]) == [(0, 1, 2), (0, 2, 3)]
    assert triangulate_polygon([(0, 0), (0, 0.5), (0, 1), (1, 1), (1, 0)]) == [(0, 2, 3), (0, 3, 4)]
