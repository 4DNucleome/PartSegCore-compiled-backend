import numpy as np
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
    split_polygon_on_repeated_edges_py,
    triangle_convex_polygon,
    triangulate_monotone_polygon_py,
    triangulate_path_edge_numpy,
    triangulate_path_edge_py,
    triangulate_polygon_numpy,
    triangulate_polygon_numpy_li,
    triangulate_polygon_py,
    triangulate_polygon_with_edge_numpy_li,
)


def test_on_segment():
    assert on_segment((0, 0), (0, 1), (0, 2))
    assert not on_segment((0, 0), (1, 1), (0, 3))


def test_orientation_basic():
    assert orientation((0, 0), (0, 1), (0, 2)) == 0
    assert orientation((0, 0), (0, 2), (0, 1)) == 0
    assert orientation((0, 2), (0, 0), (0, 1)) == 0
    assert orientation((0, 0), (0, 1), (1, 2)) == 1
    assert orientation((0, 0), (0, 1), (-1, 2)) == 2


@pytest.mark.parametrize(
    ('p1', 'p2', 'p3'),
    [
        ((0, 0), (1, 1), (2, 0)),
        ((0, 0), (1, 100), (2, 0)),
        ((0, 0), (100, 1), (1, 0)),
        ((0, 1), (1, 1), (2, 0)),
        ((0, 1), (1, 1), (2000, 0)),
    ],
)
def test_orientation_split(p1, p2, p3):
    assert orientation(p1, p2, p3) == 1
    assert orientation(p3, p2, p1) == 2


@pytest.mark.parametrize(
    ('p1', 'p2', 'p3'),
    [
        ((0, 1), (1, 0), (2, 1)),
        ((0, 1), (1, 0), (2, 0)),
        ((0, 1), (1, 0), (200, 0)),
    ],
)
def test_orientation_merge(p1, p2, p3):
    assert orientation(p1, p2, p3) == 2
    assert orientation(p3, p2, p1) == 1


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
    assert len(find_intersection_point(segment1, segment2)) == 1
    assert find_intersection_point(segment1, segment2)[0] == expected


@pytest.mark.parametrize(
    ('segment1', 'segment2', 'expected'),
    [
        (((0, 0), (2, 0)), ((1, 0), (3, 0)), [(1, 0), (2, 0)]),  # Horizontal segments
        (((1, 0), (3, 0)), ((0, 0), (2, 0)), [(1, 0), (2, 0)]),  # Horizontal segments
        (((0, 0), (0, 2)), ((0, 1), (0, 3)), [(0, 1), (0, 2)]),  # Vertical segments
        (((0, 1), (0, 3)), ((0, 0), (0, 2)), [(0, 1), (0, 2)]),  # Vertical segments
        (((0, 0), (2, 2)), ((1, 1), (3, 3)), [(1, 1), (2, 2)]),  # Diagonal segments
        (((1, 1), (3, 3)), ((0, 0), (2, 2)), [(1, 1), (2, 2)]),  # Diagonal segments
    ],
)
def test_find_intersection_points_colinear_overlap(segment1, segment2, expected):
    assert sorted(find_intersection_point(segment1, segment2)) == expected


@pytest.mark.parametrize(
    ('segment1', 'segment2', 'expected'),
    [
        (((0, 0), (3, 0)), ((1, 0), (2, 0)), [(1, 0), (2, 0)]),  # Horizontal segments
        (((1, 0), (2, 0)), ((0, 0), (3, 0)), [(1, 0), (2, 0)]),  # Horizontal segments
        (((0, 0), (0, 3)), ((0, 1), (0, 2)), [(0, 1), (0, 2)]),  # Vertical segments
        (((0, 1), (0, 2)), ((0, 0), (0, 3)), [(0, 1), (0, 2)]),  # Vertical segments
        (((0, 0), (3, 3)), ((1, 1), (2, 2)), [(1, 1), (2, 2)]),  # Diagonal segments
        (((1, 1), (2, 2)), ((0, 0), (3, 3)), [(1, 1), (2, 2)]),  # Diagonal segments
    ],
)
def test_find_intersection_point_colinear_inside(segment1, segment2, expected):
    assert sorted(find_intersection_point(segment1, segment2)) == expected


def test_is_convex():
    assert is_convex([(0, 0), (0, 1), (1, 1), (1, 0)])
    assert not is_convex([(0, 0), (0, 1), (1, 0), (1, 1)])

    assert is_convex([(0, 0), (0, 0.5), (0, 1), (1, 1), (1, 0)])
    assert not is_convex([(0, 0), (0.1, 0.5), (0, 1), (1, 1), (1, 0)])


def test_triangle_convex_polygon():
    assert triangle_convex_polygon([(0, 0), (0, 1), (1, 1), (1, 0)]) == [(0, 1, 2), (0, 2, 3)]
    assert triangle_convex_polygon([(0, 0), (0, 0.5), (0, 1), (1, 1), (1, 0)]) == [(0, 2, 3), (0, 3, 4)]


@pytest.mark.parametrize(
    ('polygon', 'expected'),
    [
        ([(0, 0), (0, 1), (1, 1), (1, 0)], [(0, 1, 2), (0, 3, 2)]),
        ([(0, 0), (0, 10), (1, 10), (1, 0)], [(0, 1, 2), (0, 3, 2)]),
        ([(0, 0), (0, 0.5), (0, 1), (1, 1), (1, 0)], [(0, 2, 3), (0, 3, 4)]),
    ],
)
def test_triangulate_polygon_py_convex(polygon, expected):
    assert is_convex(polygon)
    assert triangulate_polygon_py(polygon)[0] == expected


def _renumerate_triangles(polygon, points, triangles):
    point_num = {point: i for i, point in enumerate(polygon)}
    return [tuple(point_num[tuple(points[point])] for point in triangle) for triangle in triangles]


TEST_POLYGONS = [
    ([(0, 0), (1, 1), (0, 2), (2, 1)], [(3, 2, 1), (0, 3, 1)]),
    ([(0, 0), (0, 1), (1, 2), (2, 1), (2, 0), (1, 0.5)], [(4, 3, 5), (3, 2, 1), (5, 3, 1), (5, 1, 0)]),
    ([(0, 1), (0, 2), (1, 1.5), (2, 2), (2, 1), (1, 0.5)], [(4, 3, 2), (2, 1, 0), (4, 2, 0), (5, 4, 0)]),
    ([(0, 1), (0, 2), (1, 0.5), (2, 2), (2, 1), (1, -0.5)], [(2, 1, 0), (2, 0, 5), (4, 3, 2), (5, 4, 2)]),
    ([(0, 0), (1, 2), (2, 0), (1, 1)], [(2, 1, 3), (3, 1, 0)]),
    ([(0, 0), (0, 1), (0.5, 0.5), (1, 0), (1, 1)], [(3, 4, 2), (2, 1, 0)]),
    ([(0, 0), (1, 0), (0.5, 0.5), (0, 1), (1, 1)], [(2, 4, 3), (1, 2, 0)]),
]


@pytest.mark.parametrize(('polygon', 'expected'), TEST_POLYGONS)
def test_triangulate_polygon_py_non_convex(polygon, expected):
    triangles, points = triangulate_polygon_py(polygon)
    triangles_ = _renumerate_triangles(polygon, points, triangles)
    assert triangles_ == expected


@pytest.mark.parametrize(('polygon', 'expected'), TEST_POLYGONS)
def test_triangulate_polygon_numpy_non_convex(polygon, expected):
    triangles, points = triangulate_polygon_numpy(np.array(polygon).astype(np.float64))
    triangles_ = _renumerate_triangles(polygon, points, triangles)
    assert triangles_ == expected


@pytest.mark.parametrize(('polygon', 'expected'), TEST_POLYGONS)
def test_triangulate_polygon_numpy_li_non_convex(polygon, expected):
    triangles, points = triangulate_polygon_numpy_li([np.array(polygon).astype(np.float64)])
    triangles_ = _renumerate_triangles(polygon, points, triangles)
    assert triangles_ == expected


def test_triangulate_polygon_in_polygon_numpy():
    polygons = [
        np.array([(0, 0), (10, 0), (10, 10), (0, 10)], dtype=np.float64),
        np.array([(4, 4), (6, 4), (6, 6), (4, 6)], dtype=np.float32),
    ]
    triangles, points = triangulate_polygon_numpy_li(polygons)
    assert len(triangles) == 8
    assert len(points) == 8


def test_triangulate_polygon_segfault1():
    """Test on polygon that lead to segfault during test"""
    polygon = [
        (205.0625, 1489.83752),
        (204.212509, 1490.4751),
        (204, 1491.11255),
        (202.087509, 1493.45007),
        (201.875, 1494.7251),
        (202.300003, 1496),
        (202.300003, 1498.33752),
        (203.575012, 1499.82507),
        (204.425003, 1500.25),
        (205.0625, 1500.25),
        (205.700012, 1500.67505),
        (206.550003, 1500.67505),
        (207.1875, 1500.25),
        (208.037506, 1500.88757),
        (209.3125, 1499.82507),
        (209.525009, 1499.1875),
        (211.012512, 1497.70007),
        (210.375, 1496.42505),
        (209.525009, 1495.57507),
        (208.462509, 1495.15002),
        (208.675003, 1494.9375),
        (208.462509, 1492.8125),
        (208.037506, 1491.5376),
        (205.912506, 1489.83752),
    ]
    triangulate_polygon_py(polygon)


def test_triangulate_polygon_segfault2():
    polygon = [
        [1388.6875, 2744.4375],
        [1388.4751, 2744.6501],
        [1386.5625, 2744.6501],
        [1385.925, 2744.8625],
        [1385.5, 2745.2876],
        [1385.2876, 2747.625],
        [1385.7125, 2748.2627],
        [1385.7125, 2749.1125],
        [1386.1376, 2749.75],
        [1389.9625, 2753.7876],
        [1390.3876, 2754.6377],
        [1391.025, 2754.6377],
        [1392.0875, 2753.1501],
        [1392.3, 2753.3625],
        [1392.3, 2754.6377],
        [1392.5126, 2754.2126],
        [1392.3, 2754.0],
        [1392.3, 2751.4502],
        [1392.7251, 2750.3877],
        [1391.6626, 2748.9001],
        [1391.6626, 2747.4126],
        [1390.8125, 2745.5],
        [1390.175, 2745.2876],
        [1389.3251, 2744.4375],
    ]
    triangulate_polygon_py(polygon)


def test_triangulate_polygon_segfault3():
    polygon = [
        (1066.32507, 1794.3501),
        (1065.6875, 1794.77502),
        (1063.77502, 1794.77502),
        (1063.5625, 1794.98755),
        (1063.13757, 1794.98755),
        (1062.28748, 1795.83752),
        (1062.28748, 1797.32507),
        (1062.07507, 1797.5376),
        (1062.28748, 1797.5376),
        (1062.5, 1797.75),
        (1063.13757, 1797.75),
        (1063.5625, 1797.32507),
        (1064.625, 1797.32507),
        (1065.26257, 1797.96252),
        (1064.83752, 1797.5376),
        (1064.83752, 1796.90002),
        (1065.47498, 1796.26257),
        (1065.47498, 1796.05005),
        (1065.6875, 1795.83752),
        (1066.53748, 1795.83752),
        (1066.75, 1795.625),
        (1066.75, 1794.98755),
        (1066.96252, 1794.77502),
        (1066.75, 1794.3501),
    ]
    triangulate_polygon_py(polygon)


def test_triangulate_polygon_segfault4():
    polygon = [
        [657.6875, 2280.975],
        [657.6875, 2281.6125],
        [657.05005, 2282.25],
        [656.2, 2284.1626],
        [657.6875, 2285.8625],
        [658.11255, 2286.7126],
        [659.8125, 2288.4126],
        [659.8125, 2288.625],
        [661.9375, 2290.5376],
        [662.78754, 2290.9626],
        [664.7, 2292.2375],
        [665.3375, 2292.2375],
        [665.97504, 2291.175],
        [666.61255, 2290.5376],
        [666.61255, 2289.6875],
        [666.1875, 2288.625],
        [664.48755, 2286.925],
        [664.0625, 2286.925],
        [663.2125, 2286.5],
        [661.9375, 2284.8],
        [660.66254, 2284.8],
        [660.45, 2284.5876],
        [660.45, 2284.1626],
        [657.6875, 2281.1875],
    ]
    triangulate_polygon_py(polygon)


def test_triangulate_polygon_segfault5():
    polygon = [
        [895.05005, 2422.5],
        [894.8375, 2422.7126],
        [894.41254, 2422.7126],
        [893.98755, 2423.1375],
        [893.35004, 2423.1375],
        [892.92505, 2423.5625],
        [892.7125, 2423.5625],
        [891.4375, 2424.8376],
        [891.22504, 2424.8376],
        [892.075, 2424.8376],
        [892.5, 2425.05],
        [893.35004, 2425.9001],
        [893.775, 2426.75],
        [893.775, 2427.3875],
        [894.625, 2426.75],
        [895.6875, 2426.75],
        [896.11255, 2426.1125],
        [896.11255, 2425.05],
        [895.9, 2424.8376],
        [896.53754, 2423.7751],
        [896.53754, 2423.35],
        [896.75, 2422.925],
        [896.11255, 2422.7126],
        [895.9, 2422.5],
    ]
    triangulate_polygon_py(polygon)


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


PATH_DATA = [
    (
        [[0, 0], [0, 10], [10, 10], [10, 0]],
        True,
        False,
        10,
        [[2, 1, 0], [1, 2, 3], [4, 3, 2], [3, 4, 5], [6, 5, 4], [5, 6, 7], [8, 7, 6], [7, 8, 9]],
    ),
    (
        [[0, 0], [0, 10], [10, 10], [10, 0]],
        False,
        False,
        8,
        [[2, 1, 0], [1, 2, 3], [4, 3, 2], [3, 4, 5], [6, 5, 4], [5, 6, 7]],
    ),
    (
        [[0, 0], [0, 10], [10, 10], [10, 0]],
        True,
        True,
        14,
        [
            [2, 1, 0],
            [3, 2, 0],
            [2, 3, 4],
            [5, 4, 3],
            [6, 5, 3],
            [5, 6, 7],
            [8, 7, 6],
            [9, 8, 6],
            [8, 9, 10],
            [11, 10, 9],
            [12, 11, 9],
            [11, 12, 13],
        ],
    ),
    (
        [[0, 0], [0, 10], [10, 10], [10, 0]],
        False,
        True,
        10,
        [[2, 1, 0], [1, 2, 3], [4, 3, 2], [5, 4, 2], [4, 5, 6], [7, 6, 5], [8, 7, 5], [7, 8, 9]],
    ),
    (
        [[2, 10], [0, -5], [-2, 10], [-2, -10], [2, -10]],
        True,
        False,
        15,
        [
            [2, 1, 0],
            [1, 2, 3],
            [1, 3, 4],
            [5, 4, 3],
            [6, 5, 3],
            [5, 6, 7],
            [8, 7, 6],
            [7, 8, 9],
            [7, 9, 10],
            [11, 10, 9],
            [10, 11, 12],
            [13, 12, 11],
            [12, 13, 14],
        ],
    ),
    ([[0, 0], [0, 10]], False, False, 4, [[2, 1, 0], [1, 2, 3]]),
    ([[0, 0], [0, 10], [0, 20]], False, False, 6, [[2, 1, 0], [1, 2, 3], [4, 3, 2], [3, 4, 5]]),
    (
        [[0, 0], [0, 2], [10, 1]],
        True,
        False,
        9,
        [[2, 1, 0], [1, 2, 3], [4, 3, 2], [3, 4, 5], [6, 5, 4], [7, 6, 4], [6, 7, 8]],
    ),
    ([[0, 0], [10, 1], [9, 1.1]], False, False, 7, [[2, 1, 0], [1, 2, 3], [4, 3, 2], [3, 4, 5], [3, 5, 6]]),
    ([[9, 0.9], [10, 1], [0, 2]], False, False, 7, [[2, 1, 0], [1, 2, 3], [4, 3, 2], [3, 4, 5], [3, 5, 6]]),
    ([[0, 0], [-10, 1], [-9, 1.1]], False, False, 7, [[2, 1, 0], [1, 2, 3], [4, 3, 2], [5, 4, 2], [4, 5, 6]]),
    ([[-9, 0.9], [-10, 1], [0, 2]], False, False, 7, [[2, 1, 0], [1, 2, 3], [4, 3, 2], [5, 4, 2], [4, 5, 6]]),
]


@pytest.mark.parametrize(
    ('path', 'closed', 'bevel', 'expected', 'exp_triangles'),
    PATH_DATA,
)
@pytest.mark.parametrize('triangulate_fun', [triangulate_path_edge_py, triangulate_path_edge_numpy])
def test_triangulate_path_edge_py(path, closed, bevel, expected, exp_triangles, triangulate_fun):
    centers, offsets, triangles = triangulate_fun(np.array(path, dtype='float32'), closed=closed, bevel=bevel)
    assert centers.shape == offsets.shape
    assert centers.shape[0] == expected
    assert triangles.shape[0] == expected - 2
    triangles_li = [[int(y) for y in x] for x in triangles]
    assert triangles_li == exp_triangles


@pytest.mark.parametrize(('polygon', 'expected'), TEST_POLYGONS)
def test_triangulate_polygon_with_edge_numpy_li(polygon, expected):
    (triangles, points), (centers, offsets, _edge_triangles) = triangulate_polygon_with_edge_numpy_li(
        [np.array(polygon, dtype=np.float32)]
    )
    triangles_ = _renumerate_triangles(polygon, points, triangles)
    assert triangles_ == expected
    assert centers.shape == offsets.shape


def test_split_polygon_on_repeated_edges_py_no_split():
    res = split_polygon_on_repeated_edges_py([[0, 0], [0, 1], [1, 1], [1, 0]])
    assert len(res) == 1
    assert len(res[0]) == 4
    idx = res[0].index((0, 0))
    assert res[0][idx:] + res[0][:idx] == [(0, 0), (0, 1), (1, 1), (1, 0)]


def test_split_polygon_on_repeated_edges_py_square_in_square():
    res = split_polygon_on_repeated_edges_py(
        [[0, 0], [0, 5], [1, 5], [1, 1], [9, 1], [9, 9], [1, 9], [1, 5], [0, 5], [0, 10], [10, 10], [10, 0]]
    )
    assert len(res) == 2
    assert len(res[0]) == 5
    assert len(res[1]) == 5
    # idx = res[0].index((0, 0))
    # assert res[0][idx:] + res[0][:idx] == [(0, 0), (0, 1), (1, 1), (1, 0)]


@pytest.mark.parametrize(('split_edges', 'triangles'), [(True, 20), (False, 24)])
def test_splitting_edges(split_edges, triangles):
    polygon = np.array(
        [[0, 0], [0, 5], [1, 5], [1, 1], [9, 1], [9, 9], [1, 9], [1, 5], [0, 5], [0, 10], [10, 10], [10, 0]],
        dtype=np.float32,
    )
    triangles_ = triangulate_polygon_with_edge_numpy_li([polygon], split_edges=split_edges)[1][2]
    assert len(triangles_) == triangles


def generate_regular_polygon(n: int, reverse: bool, radius: int = 1) -> np.ndarray:
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    if reverse:
        angles = angles[::-1]
    return np.column_stack((radius * np.cos(angles), radius * np.sin(angles)))


def generate_self_intersecting_polygon(n: int, reverse: bool, radius: int = 1) -> np.ndarray:
    """Generate self-intersecting polygon with n vertices

    The polygon is generated by doubling the angle range of a regular polygon
    """
    assert n % 2 == 1, 'an odd number is required to generate a self-intersecting polygon'
    angles = np.linspace(0, 4 * np.pi, n, endpoint=False)
    if reverse:
        angles = angles[::-1]
    return np.column_stack((radius * np.cos(angles), radius * np.sin(angles)))


def rotation_matrix(angle: float) -> np.ndarray:
    """Create a 2D rotation matrix for the given angle in degrees."""
    return np.array(
        [
            [np.cos(np.radians(angle)), -np.sin(np.radians(angle))],
            [np.sin(np.radians(angle)), np.cos(np.radians(angle))],
        ]
    )


ANGLES = [0, 5, 75, 95, 355]


@pytest.mark.parametrize('angle', ANGLES, ids=str)
@pytest.mark.parametrize('n_vertex', [5, 7, 19])
@pytest.mark.parametrize('reverse', [False, True])
def test_is_convex_self_intersection(angle, n_vertex, reverse):
    p = generate_self_intersecting_polygon(n_vertex, reverse)
    rot = rotation_matrix(angle)
    data = np.dot(p, rot)
    assert not is_convex(data)


@pytest.mark.parametrize('angle', ANGLES, ids=str)
@pytest.mark.parametrize('n_vertex', [3, 4, 7, 12, 15, 20])
@pytest.mark.parametrize('reverse', [False, True])
def test_is_convex_regular_polygon(angle, n_vertex, reverse):
    poly = generate_regular_polygon(n_vertex, reverse=reverse)
    rot = rotation_matrix(angle)
    rotated_poly = np.dot(poly, rot)
    assert is_convex(rotated_poly)
