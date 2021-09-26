import numpy as np
import typing
from .euclidean_cython import calculate_euclidean as calculate_euclidean
from .fuzzy_distance import MuType as MuType, calculate_mu_array as calculate_mu_array, fuzzy_distance as fuzzy_distance
from .path_sprawl_cython import calculate_maximum as calculate_maximum, calculate_minimum as calculate_minimum
from .sprawl_utils import (
    get_closest_component as get_closest_component,
    get_maximum_component as get_maximum_component,
    get_minimum_component as get_minimum_component,
)
from typing import Any

def path_maximum_sprawl(
    data_f: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neighbourhood: np.ndarray,
    distance_cache: Any | None = ...,
    data_cache: Any | None = ...,
) -> np.ndarray: ...
def path_minimum_sprawl(
    data_f, components, components_count, neighbourhood, distance_cache: Any | None = ..., data_cache: Any | None = ...
) -> np.ndarray: ...
def euclidean_sprawl(
    data_m: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neigh_arr,
    dist_arr,
    distance_cache: Any | None = ...,
    data_cache: Any | None = ...,
) -> np.ndarray: ...
def fdt_sprawl(
    data_m: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neigh_arr,
    dist_arr,
    lower_bound,
    upper_bound,
    distance_cache: Any | None = ...,
    data_cache: Any | None = ...,
) -> np.ndarray: ...
def sprawl_component(
    data_m: np.ndarray,
    components: np.ndarray,
    component_number: int,
    neigh_arr: np.ndarray,
    dist_arr: np.ndarray,
    calculate_operator: typing.Callable,
) -> np.ndarray: ...
def distance_sprawl(
    calculate_operator,
    data_m: np.ndarray,
    components: np.ndarray,
    components_count: int,
    neigh_arr,
    dist_array,
    distance_cache: Any | None = ...,
    data_cache: Any | None = ...,
    parallel: bool = ...,
) -> np.ndarray: ...
def reverse_permutation(perm: typing.List[int]) -> typing.List[int]: ...
def relabel_with_perm(labeling: typing.List, perm: typing.List) -> typing.List: ...
def verify_cohesion(elements: typing.List[int], graph: typing.List[typing.List[int]]) -> bool: ...
def relabel_array(data: np.ndarray, perm: typing.List[int]) -> np.ndarray: ...
