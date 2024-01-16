from typing import Optional
import numpy as np

from ._napari_mapping import _zero_preserving_modulo_par, _zero_preserving_modulo_seq, _map_array_par, _map_array_seq


def _allocate_output(data: np.ndarray, num_values: int) -> np.ndarray:
    if num_values < 256:
        dtype = np.uint8
    elif num_values < 65536:
        dtype = np.uint16
    else:
        dtype = np.float32
    return np.empty_like(data, dtype=dtype)


def zero_preserving_modulo_sequential(
    data: np.ndarray, modulo_factor: int, to_zero: int, out: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Modulo plus one operation performed on values different than to_zero.
    `(n % modulo_factor) + 1 if n != to_zero else 0`

    Perform operation on each element of array sequentially.

    Parameters
    ----------
    data: np.ndarray
        data to be modified
    modulo_factor: int
        modulo factor
    to_zero: int
        value to be set to zero
    out: Optional[np.ndarray]
        output array
    """
    original_shape = data.shape
    data = data.reshape(-1)
    if out is None:
        out = _allocate_output(data, modulo_factor)
    _zero_preserving_modulo_seq(data, modulo_factor, to_zero, out)
    return out.reshape(original_shape)


def zero_preserving_modulo_parallel(
    data: np.ndarray, modulo_factor: int, to_zero: int, out: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Modulo plus one operation performed on values different than to_zero.
    `(n % modulo_factor) + 1 if n != to_zero else 0`

    Perform operation on each element of array in parallel using openmp.

    Parameters
    ----------
    data: np.ndarray
        data to be modified
    modulo_factor: int
        modulo factor
    to_zero: int
        value to be set to zero
    out: Optional[np.ndarray]
        output array
    """
    original_shape = data.shape
    data = data.reshape(-1)
    if out is None:
        out = _allocate_output(data, modulo_factor)
    _zero_preserving_modulo_par(data, modulo_factor, to_zero, out)
    return out.reshape(original_shape)


def map_array_sequential(
    data: np.ndarray, mapping: dict, default: int = 0, out: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Map values from data to values from mapping.

    Parameters
    ----------
    data: np.ndarray
        data to be modified
    mapping: dict
        dict with mapping information. Keys are values from data, values are values from output array.
        If None is provided then its value is set to default.
    default: int
        default value for not mapped values
    out: Optional[np.ndarray]
        output array
    """
    original_shape = data.shape
    data = data.reshape(-1)
    if out is None:
        out = _allocate_output(data, max(x for x in mapping.values() if x is not None) + 1)
    _map_array_seq(data, mapping, out, out.dtype.type(default))
    return out.reshape(original_shape)


def map_array_parallel(
    data: np.ndarray, mapping: dict, default: int = 0, out: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Map values from data to values from mapping.

    Parameters
    ----------
    data: np.ndarray
        data to be modified
    mapping: dict
        dict with mapping information. Keys are values from data, values are values from output array.
        If None is provided then its value is set to default.
    default: int
        default value for not mapped values
    out: Optional[np.ndarray]
        output array
    """
    original_shape = data.shape
    data = data.reshape(-1)
    if out is None:
        out = _allocate_output(data, max(x for x in mapping.values() if x is not None) + 1)
    _map_array_par(data, mapping, out, out.dtype.type(default))
    return out.reshape(original_shape)
