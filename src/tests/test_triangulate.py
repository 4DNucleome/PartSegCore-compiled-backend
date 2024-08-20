from PartSegCore_compiled_backend.triangulate import (
    on_segment,
    orientation,
    do_intersect,
    find_intersections,
    find_intersection,
    is_convex,
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
    """
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


def test_find_intersection():
    assert find_intersection(((0, 0), (2, 2)), ((0, 2), (2, 0))) == (1, 1)


def test_is_convex():
    assert is_convex([(0, 0), (0, 1), (1, 1), (1, 0)])
    assert not is_convex([(0, 0), (0, 1), (1, 0), (1, 1)])

    assert is_convex([(0, 0), (0, 0.5), (0, 1), (1, 1), (1, 0)])
    assert not is_convex([(0, 0), (0.1, 0.5), (0, 1), (1, 1), (1, 0)])


def test_triangulate_polygon():
    assert triangulate_polygon([(0, 0), (0, 1), (1, 1), (1, 0)]) == [(0, 1, 2), (0, 2, 3)]
    assert triangulate_polygon([(0, 0), (0, 10), (1, 10), (1, 0)]) == [(0, 1, 2), (0, 2, 3)]
    assert triangulate_polygon([(0, 0), (0, 0.5), (0, 1), (1, 1), (1, 0)]) == [(0, 2, 3), (0, 3, 4)]
