import typing
from enum import Enum

import numpy as np

class MuType(Enum):
    base_mu = 1
    reflection_mu = 2
    two_object_mu = 3

def calculate_mu(
    image: np.ndarray,
    lower_bound: float,
    upper_bound: float,
    type_: MuType,
    mask: typing.Optional[np.ndarray] = None,
    lower_mid_bound: float = 0,
    upper_mid_bound: float = 0,
): ...

class PyMSO:
    def set_image(
        self,
        image: np.ndarray,
        image_types: float,
        lower_bound: float,
        upper_bound: float,
        type_: MuType,
        mask: typing.Optional[np.ndarray] = None,
        lower_mid_bound: float = 0,
        upper_mid_bound: float = 0,
    ): ...
    def set_mu_array(self, mu: np.ndarray): ...
    def set_components(self, components: np.ndarray, component_num: typing.Optional[int] = None): ...
    def set_neighbourhood(self, neighbourhood: np.ndarray, distances: np.ndarray): ...
    def calculate_FDT(self) -> np.ndarray: ...
    def optimum_erosion_calculate(
        self, fdt_array: np.ndarray, components_arr: np.ndarray, sprawl_area: np.ndarray
    ) -> np.ndarray: ...
    def constrained_dilation(
        self, fdt_array: np.ndarray, components_arr: np.ndarray, sprawl_area: np.ndarray
    ) -> np.ndarray: ...
    def run_MSO(self, step_limits: int = 1, count_steps_factor: int = 3): ...
    def steps_done(self) -> int: ...
    def set_components_num(self, num: int): ...
    def get_result_catted(self) -> np.ndarray: ...
    def get_fdt(self) -> np.ndarray: ...
    def set_use_background(self, use: bool): ...

def calculate_mu_mid(image: np.ndarray, lower_bound: float, mid_point: float, upper_bound: float): ...
