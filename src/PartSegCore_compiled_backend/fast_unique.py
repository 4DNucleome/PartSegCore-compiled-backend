import numpy as np

from ._fast_unique import unique2d, unique3d


def label_unique(array, numpy_fallback=False):
    try:
        if array.ndim == 2:
            return unique2d(array)
        elif array.ndim == 3:
            return unique3d(array)
        elif numpy_fallback:
            return np.unique(array)
        raise RuntimeError('Array must be 2d or 3d')
    except TypeError:
        if numpy_fallback:
            return np.unique(array)
        else:
            raise
