from typing import Optional
import numpy as np

def get_maximum_component(
    components: np.ndarray,
    data_mask: np.ndarray,
    paths: np.ndarray,
    components_translation: np.ndarray,
    num_of_components: Optional[int] = None,
) -> np.ndarray: ...
def get_minimum_component(
    components: np.ndarray,
    data_mask: np.ndarray,
    paths: np.ndarray,
    components_translation: np.ndarray,
    num_of_components: Optional[int] = None,
) -> np.ndarray: ...
def get_closest_component(
    components: np.ndarray,
    data_mask: np.ndarray,
    distances: np.ndarray,
    components_translation: np.ndarray,
    num_of_components: Optional[int] = None,
) -> np.ndarray: ...
