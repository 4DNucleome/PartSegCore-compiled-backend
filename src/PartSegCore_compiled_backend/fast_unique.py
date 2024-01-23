"""
Module with fast unique calculation.
This implementation is faster than numpy.unique and comparable to
pandas.unique(data.ravel())` for data with continuous memory.
If data are not continuous memory this function is faster than pandas.unique.

For benchmark please see https://github.com/4DNucleome/PartSegCore-compiled-backend/blob/master/notebooks/performance_label_unique.ipynb
"""
import numpy as np

from ._fast_unique import unique1d, unique2d, unique3d


def label_unique(array, numpy_fallback=False) -> np.ndarray:
    """
    Calculate unique values in array.

    Parameters
    ----------
    array : np.ndarray
        array to calculate unique values
    numpy_fallback : bool
        if True allow using numpy.unique if cython version is not available
        otherwise raise RuntimeError

    Returns
    -------
    np.ndarray
        array of unique values
    """
    if array.ndim == 1:
        return unique1d(array)
    if array.ndim == 2:
        return unique2d(array)
    if array.ndim == 3:
        return unique3d(array)
    if numpy_fallback:
        return np.unique(array)
    raise RuntimeError('Array must be 1d, 2d or 3d')
