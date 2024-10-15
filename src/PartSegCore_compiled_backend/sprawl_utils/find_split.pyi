import typing
from typing import Any

import numpy as np

from PartSegCore_compiled_backend.sprawl_utils.euclidean_cython import calculate_euclidean as calculate_euclidean
from PartSegCore_compiled_backend.sprawl_utils.fuzzy_distance import MuType as MuType
from PartSegCore_compiled_backend.sprawl_utils.fuzzy_distance import calculate_mu_array as calculate_mu_array
from PartSegCore_compiled_backend.sprawl_utils.fuzzy_distance import fuzzy_distance as fuzzy_distance
from PartSegCore_compiled_backend.sprawl_utils.path_sprawl_cython import calculate_maximum as calculate_maximum
from PartSegCore_compiled_backend.sprawl_utils.path_sprawl_cython import calculate_minimum as calculate_minimum
from PartSegCore_compiled_backend.sprawl_utils.sprawl_utils import get_closest_component as get_closest_component
from PartSegCore_compiled_backend.sprawl_utils.sprawl_utils import get_maximum_component as get_maximum_component
from PartSegCore_compiled_backend.sprawl_utils.sprawl_utils import get_minimum_component as get_minimum_component

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
def reverse_permutation(perm: list[int]) -> list[int]: ...
def relabel_with_perm(labeling: list, perm: list) -> list: ...
def verify_cohesion(elements: list[int], graph: list[list[int]]) -> bool: ...
def relabel_array(data: np.ndarray, perm: list[int]) -> np.ndarray: ...
