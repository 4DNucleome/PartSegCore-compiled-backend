import numpy as np
import numpy.testing as npt
import pytest
from PartSegCore_compiled_backend.napari_mapping import (
    zero_preserving_modulo_parallel,
    zero_preserving_modulo_sequential,
    map_array_parallel,
    map_array_sequential,
)


def _zero_preserving_modulo_numpy(values: np.ndarray, n: int, dtype: np.dtype, to_zero: int = 0) -> np.ndarray:
    res = ((values - 1) % n + 1).astype(dtype)
    res[values == to_zero] = 0
    return res


def _map_array(values: np.ndarray, map_dict: dict, dtype: np.dtype) -> np.ndarray:
    return np.array([map_dict.get(x, 0) for x in values.flatten()], dtype=dtype).reshape(values.shape)


DATA_LI = [
    np.arange(100, dtype=np.uint32),
    np.arange(100, dtype=np.int32),
    np.arange(100, dtype=np.uint64),
    np.arange(100, dtype=np.int64),
    np.arange(100, dtype=np.uint32).reshape((10, 10)),
]
DATA_IDS = ['uint32', 'int32', 'uint64', 'int64', 'uint32_2d']


@pytest.mark.parametrize('func', [zero_preserving_modulo_parallel, zero_preserving_modulo_sequential])
@pytest.mark.parametrize('data', DATA_LI, ids=DATA_IDS)
def test_zero_preserving_modulo(func, data):
    out = func(data, 10, 0)
    assert out.dtype == np.uint8
    npt.assert_array_equal(out, _zero_preserving_modulo_numpy(data, 10, np.uint8, 0))


@pytest.mark.parametrize('func', [map_array_parallel, map_array_sequential])
@pytest.mark.parametrize('data', DATA_LI, ids=DATA_IDS)
def test_map_array(func, data):
    map_dkt = {1: 1, 2: 2, 3: 1, 4: 3, 5: 2}
    out = func(data, map_dkt, 0)
    assert out.dtype == np.uint8
    npt.assert_array_equal(out, _map_array(data, map_dkt, np.uint8))
