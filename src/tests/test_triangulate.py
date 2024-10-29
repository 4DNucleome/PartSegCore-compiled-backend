import pytest

from PartSegCore_compiled_backend.triangulate import (
    do_intersect,
    find_intersection_point,
    find_intersection_points_py,
    find_intersections,
    is_convex,
    on_segment,
    orientation,
    segment_left_to_right_comparator,
    triangle_convex_polygon,
    triangulate_monotone_polygon_py,
    triangulate_polygon,
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
    ('segments', 'expected'),
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
    ('segment1', 'segment2', 'expected'),
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


@pytest.mark.parametrize(
    ('segment1', 'segment2'),
    [
        # Touching segments by top point
        (((0, 0), (1, 1)), ((1, 1), (2, 0))),
        # Touching segments by bottom point
        (((0, 1), (1, 0)), ((1, 0), (2, 1))),
        # Touching segments by bottom point with different length, longer left
        (((-1, 2), (1, 0)), ((1, 0), (2, 1))),
        # Touching segments by bottom point with different length, longer right
        (((0, 1), (1, 0)), ((1, 0), (3, 2))),
        # Parallel vertical segments
        (((0, 0), (0, 1)), ((1, 0), (1, 1))),
        # Parallel horizontal segments different length
        (((0, 0), (0, 2)), ((1, -1), (1, 1))),
        # horizontal segments
        (((0, 0), (1, 0)), ((2, 0), (3, 0))),
        # Two segments with top on same line
        (((0, 0), (1, 1)), ((3, 1), (2, -1))),
        # One horizontal segment on right and top
        (((0, 0), (1, 1)), ((2, 1), (3, 1))),
        # One horizontal segment on right and middle
        (((0, 0), (1, 1)), ((2, 0.5), (3, 0.5))),
        # One horizontal segment on right and bottom
        (((0, 0), (1, 1)), ((2, 0), (3, 0))),
        # One horizontal segment on left and top
        (((0, 1), (1, 1)), ((2, 0), (3, 2))),
        # One horizontal segment on left and middle
        (((0, 0.5), (1, 0.5)), ((2, 0), (3, 2))),
        # One horizontal segment on left and bottom
        (((0, 0), (1, 0)), ((2, 0), (3, 2))),
        # left horizontal, right oblique, bottom merge
        (((0, 0), (1, 0)), ((1, 0), (2, 2))),
        # left oblique, right horizontal, bottom merge
        (((0, 1), (1, 0)), ((1, 0), (2, 0))),
        # left horizontal, right oblique, top merge
        (((0, 1), (1, 1)), ((1, 1), (2, 0))),
        # left oblique, right horizontal, top merge
        (((0, 0), (1, 1)), ((1, 1), (2, 1))),
    ],
    ids=[
        'Touching segments by top point',
        'Touching segments by bottom point',
        'Touching segments by bottom point with different length, longer left',
        'Touching segments by bottom point with different length, longer right',
        'Parallel vertical segments',
        'Parallel horizontal segments different length',
        'horizontal segments',
        'Two segments with top on same line',
        'One horizontal segment on right and top',
        'One horizontal segment on right and middle',
        'One horizontal segment on right and bottom',
        'One horizontal segment on left and top',
        'One horizontal segment on left and middle',
        'One horizontal segment on left and bottom',
        'left horizontal, right oblique, bottom merge',
        'left oblique, right horizontal, bottom merge',
        'left horizontal, right oblique, top merge',
        'left oblique, right horizontal, top merge',
    ],
)
def test_segment_left_to_right_comparator(segment1, segment2):
    assert segment_left_to_right_comparator(segment1, segment2)
    assert not segment_left_to_right_comparator(segment2, segment1)


def test_find_intersection_points_py_cross():
    r"""
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
    assert find_intersection_points_py([(0, 0), (1, 1), (1, 0), (0, 1)]) == [
        (0, 0),
        (0.5, 0.5),
        (1, 1),
        (1, 0),
        (0.5, 0.5),
        (0, 1),
    ]


def test_find_intersection_points_py_cross_intersect_in_point():
    r"""
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
    assert find_intersection_points_py([(0, 0), (0.5, 0.5), (1, 1), (1, 0), (0, 1)]) == [
        (0, 0),
        (0.5, 0.5),
        (1, 1),
        (1, 0),
        (0.5, 0.5),
        (0, 1),
    ]


@pytest.mark.parametrize(
    ('polygon', 'expected'),
    [
        (
            ((1, 2), (1, 0), [(0, 1)], [(2, 1)]),
            [[(2.0, 1.0), (1.0, 2.0), (0.0, 1.0)], [(1.0, 0.0), (2.0, 1.0), (0.0, 1.0)]],
        ),
        (
            ((5, 2), (5, 0), [(0, 1)], [(2, 1)]),
            [[(2.0, 1.0), (5.0, 2.0), (0.0, 1.0)], [(5.0, 0.0), (2.0, 1.0), (0.0, 1.0)]],
        ),
        (
            ((1, 3), (1, 0), [(0, 2), (0, 1)], [(2, 2), (2, 1)]),
            [
                [(2.0, 2.0), (1.0, 3.0), (0.0, 2.0)],
                [(2.0, 1.0), (2.0, 2.0), (0.0, 2.0)],
                [(2.0, 1.0), (0.0, 2.0), (0.0, 1.0)],
                [(1.0, 0.0), (2.0, 1.0), (0.0, 1.0)],
            ],
        ),
        (
            ((0, 4), (0, 0), [], [(2, 3), (3, 2), (2, 1)]),
            [
                [(3.0, 2.0), (2.0, 3.0), (0.0, 4.0)],
                [(2.0, 1.0), (3.0, 2.0), (0.0, 4.0)],
                [(2.0, 1.0), (0.0, 4.0), (0.0, 0.0)],
            ],
        ),
        (
            ((0, 4), (0, 0), [], [(2, 3), (1, 2), (2, 1)]),
            [
                [(1.0, 2.0), (2.0, 3.0), (0.0, 4.0)],
                [(1.0, 2.0), (0.0, 4.0), (0.0, 0.0)],
                [(2.0, 1.0), (1.0, 2.0), (0.0, 0.0)],
            ],
        ),
    ],
)
def test_triangulate_monotone_polygon_py(polygon, expected):
    assert triangulate_monotone_polygon_py(*polygon) == expected
