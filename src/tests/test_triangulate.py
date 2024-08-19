from PartSegCore_compiled_backend.triangulate import on_segment, orientation, do_intersect, find_intersections


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
