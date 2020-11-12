from itertools import product

import numpy as np
import pytest

from PartSegCore_compiled_backend.color_image_cython import calculate_borders, calculate_borders2d


@pytest.mark.parametrize("label_type", [np.uint8, np.uint16, np.uint32])
def test_calculate_borders(label_type):
    layers = np.zeros((1, 10, 10, 10), dtype=label_type)
    layers[:, 3:-3, 3:-3, 3:-3] = 1
    res = calculate_borders(layers, 0, False)
    expected = layers.copy()
    expected[:, 4:-4, 4:-4, 4:-4] = 0
    assert np.all(res == expected)
    res = calculate_borders(layers, 0, True)
    expected = layers.copy()
    expected[:, :, 4:-4, 4:-4] = 0
    assert np.all(res == expected)
    res = calculate_borders2d(layers, 0)
    assert np.all(res == expected)


@pytest.mark.parametrize("label_type", [np.uint8, np.uint16, np.uint32])
def test_calculate_borders_thick(label_type):
    layers = np.zeros((1, 16, 16, 16), dtype=label_type)
    layers[:, 3:-3, 3:-3, 3:-3] = 1
    res = calculate_borders(layers, 1, False)
    expected = np.zeros((1, 16, 16, 16), dtype=np.uint8)
    expected[:, 2:-2, 2:-2, 2:-2] = 1
    expected[:, 5:-5, 5:-5, 5:-5] = 0
    for c1, c2 in product([2, -3], repeat=2):
        for x in range(3):
            cord = [c1, c2]
            cord.insert(x, slice(None))
            expected[tuple([0] + cord)] = 0

    assert np.all(res == expected)


@pytest.mark.parametrize("label_type", [np.uint8, np.uint16, np.uint32])
def test_calculate_borders_thick2d(label_type):
    layers = np.zeros((1, 16, 16, 16), dtype=label_type)
    layers[:, 3:-3, 3:-3, 3:-3] = 1
    res1 = calculate_borders(layers, 1, True)
    res2 = calculate_borders2d(layers, 1)
    assert np.all(res1 == res2)
    expected = np.zeros((1, 16, 16, 16), dtype=np.uint8)
    expected[:, 3:-3, 2:-2, 2:-2] = 1
    expected[:, :, 5:-5, 5:-5] = 0
    for c1, c2 in product([2, -3], repeat=2):
        expected[0, :, c1, c2] = 0

    assert np.all(res1 == expected)
    assert np.all(res2 == expected)
