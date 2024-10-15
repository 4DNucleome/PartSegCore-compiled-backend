import numpy as np

def calculate_euclidean(
    object_area: np.ndarray,
    base_object: np.ndarray,
    neighbourhood: np.ndarray,
    distance: np.ndarray,
    distance_cache: np.ndarray | None = None,
) -> np.ndarray: ...
