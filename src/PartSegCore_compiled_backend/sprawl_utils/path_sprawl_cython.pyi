from typing import Optional
import numpy as np

def calculate_maximum(
    object_area: np.ndarray, base_object: np.ndarray, neighbourhood: np.ndarray, result: Optional[np.ndarray] = None
) -> np.ndarray: ...
def calculate_minimum(
    object_area: np.ndarray, base_object: np.ndarray, neighbourhood: np.ndarray, result: Optional[np.ndarray] = None
) -> np.ndarray: ...
