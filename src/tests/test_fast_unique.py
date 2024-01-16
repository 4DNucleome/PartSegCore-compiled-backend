import numpy as np
import pytest

from PartSegCore_compiled_backend.fast_unique import label_unique


@pytest.mark.parametrize(
    'array',
    [
        np.random.default_rng().integers(0, 100, size=(1000,)),
        np.random.default_rng().integers(0, 100, size=(100, 100)),
        np.random.default_rng().integers(0, 100, size=(10, 100, 100)),
    ],
    ids=['1d', '2d', '3d'],
)
def test_label_unique(array):
    assert np.array_equal(np.unique(array), label_unique(array))


def test_label_unique_fallback():
    data = np.random.default_rng().integers(0, 100, size=(10, 10, 10, 10))
    with pytest.raises(RuntimeError, match='Array must be'):
        label_unique(data)
    assert np.array_equal(np.unique(data), label_unique(data, numpy_fallback=True))
