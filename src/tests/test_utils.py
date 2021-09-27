import numpy as np
import pytest

from PartSegCore_compiled_backend.utils import calc_bounds


@pytest.mark.parametrize("ndim", range(2, 6))
@pytest.mark.parametrize("dtype", [np.uint8, np.uint16, np.uint32])
def test_calc_bounds(ndim, dtype):
    data = np.zeros((20,) + (10,) * (ndim - 1), dtype=dtype)
    slice_arr = [slice(2, 8) for _ in range(ndim)]
    data[tuple(slice_arr)] = 1
    slice_arr[0] = slice(12, 18)
    data[tuple(slice_arr)] = 3
    low_bound, upper_bound = calc_bounds(data)
    assert low_bound.shape == (4, ndim)
    assert np.all(low_bound[0] == 25)
    assert np.all(upper_bound[0] == -1)
    assert np.all(low_bound[1] == 2)
    assert np.all(upper_bound[1] == 7)
    assert np.all(low_bound[2] == 25)
    assert np.all(upper_bound[2] == -1)
    assert np.all(low_bound[3][1:] == 2)
    assert np.all(upper_bound[3][1:] == 7)
    assert low_bound[3][0] == 12
    assert upper_bound[3][0] == 17
