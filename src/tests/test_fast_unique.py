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
    with pytest.raises(RuntimeError, match='Array must be'):
        label_unique(np.random.default_rng().integers(0, 100, size=(10, 10, 10, 10)))
    assert np.array_equal(
        np.unique(np.random.default_rng().integers(0, 100, size=(10, 10, 10, 10))),
        label_unique(np.random.default_rng().integers(0, 100, size=(10, 10, 10, 10)), True),
    )
