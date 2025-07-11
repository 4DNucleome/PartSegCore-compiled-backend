"""
This module contains function to perform multiple component sprawl (watershed like)
"""

from __future__ import annotations

import os
import typing
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial

import numpy as np

from PartSegCore_compiled_backend.sprawl_utils.euclidean_cython import calculate_euclidean
from PartSegCore_compiled_backend.sprawl_utils.fuzzy_distance import MuType, calculate_mu_array, fuzzy_distance
from PartSegCore_compiled_backend.sprawl_utils.path_sprawl_cython import calculate_maximum, calculate_minimum
from PartSegCore_compiled_backend.sprawl_utils.sprawl_utils import (
    get_closest_component,
    get_maximum_component,
    get_minimum_component,
)


def _allocate_cache(distance_cache: np.ndarray | None, data_shape: tuple):
    if distance_cache is None or distance_cache.shape[0] < 3:
        cache_size = 5
        distance_cache = np.zeros((cache_size + 1, *data_shape), dtype=np.float64)
    else:
        cache_size = distance_cache.shape[0] - 1
    components_numbers_translate = np.zeros(cache_size + 1, dtype=np.uint32)
    return distance_cache, cache_size, components_numbers_translate


def path_maximum_sprawl(
    data_f: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neighbourhood: np.ndarray,
    distance_cache: np.ndarray | None = None,
    data_cache: np.ndarray | None = None,
) -> np.ndarray:
    """
    Calculate sprawl in respect to brightness. Distance between voxels is the maximum brightness on
    all paths connecting them.

    Parameters
    ----------
    data_f : np.ndarray
        Array with brightness.
    components : np.ndarray
        Core components as labeled array.
    components_count : int
        Number of components.
    neighbourhood : np.ndarray
        Information about neighbourhood, shift array.
    distance_cache : np.ndarray, optional
        Cache array for reducing memory allocation. Shape of data_m, dtype is np.float64.
    data_cache : np.ndarray, optional
        Cache array for reducing memory allocation. Shape of data_m, dtype is data_m.dtype
        (use np.uint8 instead of np.bool)

    Returns
    -------
    np.ndarray
        Array with updated labels.
    """

    if data_cache is None:
        data_cache = np.zeros(data_f.shape, data_f.dtype)
    if components_count == 1:
        np.copyto(data_cache, data_f)
        tmp = calculate_maximum(data_cache, (components == 1).astype(np.uint8), neighbourhood)
        components[tmp > 0] = 1
        return components
    base_components = np.copy(components)
    masked_area = ((data_f > 0) * (base_components == 0)).astype(np.uint8)
    distance_cache, cache_size, components_numbers_translate = _allocate_cache(distance_cache, data_f.shape)
    distance_cache[:] = -np.inf
    current_in_cache = 0
    for component in range(1, components_count + 1):
        np.copyto(data_cache, data_f)
        data_cache[(base_components > 0) * (base_components != component)] = 0
        current_in_cache += 1
        distance_cache[current_in_cache] = calculate_maximum(
            data_cache, (components == component).astype(np.uint8), neighbourhood
        )
        components_numbers_translate[current_in_cache] = component
        if current_in_cache == cache_size:
            components = get_maximum_component(
                components, masked_area, distance_cache, components_numbers_translate, current_in_cache
            )
            current_in_cache = 0
    if current_in_cache > 0:
        components = get_maximum_component(
            components, masked_area, distance_cache, components_numbers_translate, current_in_cache
        )
    return components


def path_minimum_sprawl(
    data_f: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neighbourhood: np.ndarray,
    distance_cache: np.ndarray | None = None,
    data_cache: np.ndarray | None = None,
) -> np.ndarray:
    """
    Calculate sprawl in respect to brightness. Distance between voxels is the maximum brightness on
    all paths connecting them.

    Parameters
    ----------
    data_f : np.ndarray
        Array with brightness.
    components : np.ndarray
        Core components as labeled array.
    components_count : int
        Number of components.
    neighbourhood : np.ndarray
        Information about neighbourhood, shift array.
    distance_cache : np.ndarray, optional
        Cache array for reducing memory allocation. Shape of data_m, dtype is np.float64.
    data_cache : np.ndarray, optional
        Cache array for reducing memory allocation. Shape of data_m, dtype is data_m.dtype (use np.uint8 instead of np.bool).

    Returns
    -------
    np.ndarray
        Array with updated labels.
    """
    maximum = data_f.max()
    if data_cache is None:
        data_cache = np.zeros(data_f.shape, data_f.dtype)
    if components_count == 1:
        np.copyto(data_cache, data_f)
        tmp = calculate_minimum(data_cache, (components == 1).astype(np.uint8), neighbourhood)
        components[tmp < maximum] = 1
        return components
    base_components = np.copy(components)
    masked_area = ((data_f > 0) * (base_components == 0)).astype(np.uint8)
    distance_cache, cache_size, components_numbers_translate = _allocate_cache(distance_cache, data_f.shape)
    distance_cache[:] = np.inf
    current_in_cache = 0
    for component in range(1, components_count + 1):
        np.copyto(data_cache, data_f)
        data_cache[(base_components > 0) * (base_components != component)] = 0
        current_in_cache += 1
        distance_cache[current_in_cache] = calculate_minimum(
            data_cache, (components == component).astype(np.uint8), neighbourhood
        )
        components_numbers_translate[current_in_cache] = component
        if current_in_cache == cache_size:
            components = get_minimum_component(
                components, masked_area, distance_cache, components_numbers_translate, current_in_cache
            )
            current_in_cache = 0
    if current_in_cache > 0:
        components = get_minimum_component(
            components, masked_area, distance_cache, components_numbers_translate, current_in_cache
        )
    return components


def euclidean_sprawl(
    data_m: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neigh_arr,
    dist_arr,
    distance_cache=None,
    data_cache=None,
) -> np.ndarray:
    """
    Calculate euclidean sprawl (watershed)

    Parameters
    ----------
    data_m : np.ndarray
        Area for sprawl.
    components : np.ndarray
        Core components as labeled array.
    components_count : int
        Number of components.
    neigh_arr : object
        Information about neighbourhood, shift array.
    dist_arr : object
        Information about neighbourhood, distance array for shifted position.
    distance_cache : np.ndarray, optional
        Cache array for reducing memory allocation. Shape of data_m, dtype is np.float64.
    data_cache : np.ndarray, optional
        Cache array for reducing memory allocation. Shape of data_m, dtype is data_m.dtype
        (use np.uint8 instead of np.bool).

    Returns
    -------
    np.ndarray
        Array with updated labels.
    """
    #  return calculate_euclidean_iterative(data_m, components, neigh_arr, dist_arr, components_count)
    return distance_sprawl(
        calculate_euclidean, data_m, components, components_count, neigh_arr, dist_arr, distance_cache, data_cache
    )


def fdt_sprawl(
    data_m: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neigh_arr,
    dist_arr,
    lower_bound,
    upper_bound,
    distance_cache=None,
    data_cache=None,
) -> np.ndarray:
    """
    Calculate the FDT (Fuzzy Distance Transform) sprawl.

    Parameters
    ----------
    data_m : np.ndarray
        The sprawl area.
    components : np.ndarray
        Core components as labeled array.
    components_count : int
        Number of components.
    neigh_arr : array-like
        Information about the neighborhood, shift array.
    dist_arr : array-like
        Information about the neighborhood, distance array for shifted position.
    lower_bound : float
        Lower bound for calculating the mu value.
    upper_bound : float
        Upper bound for calculating the mu value.
    distance_cache : np.ndarray, optional
        Cache array for reducing memory allocation. Shape of data_m, dtype is np.float64.
    data_cache : np.ndarray, optional
        Cache array for reducing memory allocation. Shape of data_m, dtype is data_m.dtype
        (use np.uint8 instead of np.bool).

    Returns
    -------
    np.ndarray
        Array with updated labels.
    """
    if lower_bound > upper_bound:
        mu_array = 1 - calculate_mu_array(data_m, upper_bound, lower_bound, MuType.reflection_mu)
    else:
        mu_array = calculate_mu_array(data_m, lower_bound, upper_bound, MuType.reflection_mu)
    return distance_sprawl(
        partial(fuzzy_distance, mu_array=mu_array),
        data_m,
        components,
        components_count,
        neigh_arr,
        dist_arr,
        distance_cache,
        data_cache,
    )


def sprawl_component(
    data_m: np.ndarray,
    components: np.ndarray,
    component_number: int,
    neigh_arr: np.ndarray,
    dist_arr: np.ndarray,
    calculate_operator: typing.Callable,
) -> np.ndarray:
    """
    Calculate the sprawl for a single component.

    Parameters
    ----------
    data_m : np.ndarray
        The sprawl area.
    components : np.ndarray
        Array of components.
    component_number : int
        Chosen component number.
    neigh_arr : np.ndarray
        Information about the neighborhood, shift array.
    dist_arr : np.ndarray
        Information about the neighborhood, distance array for shifted position.
    calculate_operator : typing.Callable
        Function to be called to calculate sprawl.

    Returns
    -------
    np.ndarray
        Array with updated labels.
    """
    data_cache = np.copy(data_m)
    data_cache[(components > 0) * (components != component_number)] = 0
    return calculate_operator(
        data_cache, np.array(components == component_number).astype(np.uint8), neigh_arr, dist_arr
    )


def distance_sprawl(
    calculate_operator,
    data_m: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neigh_arr,
    dist_array,
    distance_cache=None,
    data_cache=None,
    parallel=False,
) -> np.ndarray:
    if data_m.dtype == bool:
        data_m = data_m.astype(np.uint8)
    if data_cache is None:
        data_cache = np.zeros(data_m.shape, data_m.dtype)
    if components_count == 1:
        np.copyto(data_cache, data_m)
        tmp = calculate_operator(data_cache, (components == 1).astype(np.uint8), neigh_arr, dist_array)
        components[tmp < 2**17] = 1
        return components
    base_components = np.copy(components)
    masked_area = ((data_m > 0) * (base_components == 0)).astype(np.uint8)
    distance_cache, cache_size, components_numbers_translate = _allocate_cache(distance_cache, data_m.shape)
    distance_cache[:] = np.inf
    current_in_cache = 0
    if parallel:
        workers_num = os.cpu_count()
        if not workers_num:
            workers_num = 4
        with ThreadPoolExecutor(max_workers=workers_num) as executor:
            result_info = {
                executor.submit(
                    sprawl_component,
                    data_m,
                    base_components,
                    component_number,
                    neigh_arr,
                    dist_array,
                    calculate_operator,
                ): component_number
                for component_number in range(1, components_count + 1)
            }
            for result in as_completed(result_info):
                component_number = result_info[result]
                current_in_cache += 1
                components_numbers_translate[current_in_cache] = component_number
                distance_cache[current_in_cache] = result.result()
                if current_in_cache == cache_size:
                    components = get_closest_component(
                        components, masked_area, distance_cache, components_numbers_translate, current_in_cache
                    )
                    current_in_cache = 0

    else:
        for component in range(1, components_count + 1):
            np.copyto(data_cache, data_m)
            data_cache[(base_components > 0) * (base_components != component)] = 0
            current_in_cache += 1
            distance_cache[current_in_cache] = calculate_operator(
                data_cache,
                (components == component).astype(np.uint8),
                neigh_arr,
                dist_array,
                distance_cache=distance_cache[0],
            )
            components_numbers_translate[current_in_cache] = component
            if current_in_cache == cache_size:
                components = get_closest_component(
                    components, masked_area, distance_cache, components_numbers_translate, current_in_cache
                )
                current_in_cache = 0
    if current_in_cache > 0:
        components = get_closest_component(
            components, masked_area, distance_cache, components_numbers_translate, current_in_cache
        )
    return components


def reverse_permutation(perm: list[int]) -> list[int]:
    rev = [0] * len(perm)
    for i, x in enumerate(perm, 1):
        rev[x - 1] = i
    return rev


def relabel_with_perm(labeling: list, perm: list) -> list:
    perm = reverse_permutation(perm)
    return [perm[x] for x in labeling]


def verify_cohesion(elements: list[int], graph: list[list[int]]) -> bool:
    start = elements[0]
    elements_set = set(elements)
    elements_set.remove(start)
    queue = list(graph[start])
    while queue and elements_set:
        el = queue.pop()
        if el in elements_set:
            elements_set.remove(el)
            queue.extend(graph[el])
    return len(elements_set) == 0


def relabel_array(data: np.ndarray, perm: list[int]) -> np.ndarray:
    result = np.zeros(data.shape, dtype=data.dtype)
    for i, val in enumerate(perm):
        result[data == i + 1] = val
    return result
