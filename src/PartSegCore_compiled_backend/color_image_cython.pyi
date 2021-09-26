import numpy as np

def calculate_borders(labels: np.ndarray, border_thick: int, per_layer: bool = True) -> np.ndarray: ...
def calculate_borders2d(labels: np.ndarray, border_thick: int) -> np.ndarray: ...
def color_grayscale(
    cmap: np.ndarray, image: np.ndarray, min_val: float, max_val: float, single_channel: int = -1
) -> np.ndarray: ...
def add_labels(
    image: np.ndarray,
    labels: np.ndarray,
    overlay: float,
    only_border: int,
    border_thick: int,
    use_labels: np.ndarray,
    label_colors: np.ndarray,
) -> np.ndarray: ...
