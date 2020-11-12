from itertools import product

import numpy as np
import pytest

from PartSegCore_compiled_backend.sprawl_utils.euclidean_cython import calculate_euclidean


@pytest.fixture
def cube_data():
    data = np.zeros((10, 10, 10), dtype=np.uint8)
    data[3:-3, 3:-3, 3:-3] = 1
    return data


@pytest.fixture
def neigh():
    return np.array([x for x in product(range(-1, 2), repeat=3) if x != (0, 0, 0)], dtype=np.int8)


def test_calculate_euclidean(cube_data, neigh):
    res = calculate_euclidean(
        np.ones(cube_data.shape, dtype=np.uint8), cube_data, neigh, np.max(np.abs(neigh), axis=1).astype(np.float64)
    )
    expected = np.ones(cube_data.shape, dtype=np.uint8) * 3
    expected[1:-1, 1:-1, 1:-1] = 2
    expected[2:-2, 2:-2, 2:-2] = 1
    expected[3:-3, 3:-3, 3:-3] = 0
    assert np.all(res == expected)


def test_calculate_euclidean_mask(cube_data, neigh):
    mask = np.zeros(cube_data.shape, dtype=np.uint8)
    mask[1:-1, 1:-1, 1:-1] = 1
    res = calculate_euclidean(mask, cube_data, neigh, np.max(np.abs(neigh), axis=1).astype(np.float64))
    expected = np.zeros(cube_data.shape, dtype=np.float64)
    expected[:] = np.inf
    expected[1:-1, 1:-1, 1:-1] = 2
    expected[2:-2, 2:-2, 2:-2] = 1
    expected[3:-3, 3:-3, 3:-3] = 0
    assert np.all(res == expected)
