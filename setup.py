import os

import numpy as np
from setuptools import Extension, setup

current_dir = os.path.dirname(os.path.abspath(__file__))
package_dir = os.path.join(current_dir, "src")

cpp_standard = ["-std=c++11"]
sprawl_utils_path = [os.path.join(package_dir, "PartSegCore_compiled_backend", "sprawl_utils")]


extensions = [
    Extension(
        "PartSegCore_compiled_backend.sprawl_utils.euclidean_cython",
        sources=["src/PartSegCore_compiled_backend/sprawl_utils/euclidean_cython.pyx"],
        include_dirs=[np.get_include()] + sprawl_utils_path,
        language="c++",
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
    ),
    Extension(
        "PartSegCore_compiled_backend.sprawl_utils.path_sprawl_cython",
        sources=["src/PartSegCore_compiled_backend/sprawl_utils/path_sprawl_cython.pyx"],
        include_dirs=[np.get_include()] + sprawl_utils_path,
        language="c++",
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
    ),
    Extension(
        "PartSegCore_compiled_backend.sprawl_utils.sprawl_utils",
        sources=["src/PartSegCore_compiled_backend/sprawl_utils/sprawl_utils.pyx"],
        include_dirs=[np.get_include()] + sprawl_utils_path,
        language="c++",
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
    ),
    Extension(
        "PartSegCore_compiled_backend.sprawl_utils.fuzzy_distance",
        sources=["src/PartSegCore_compiled_backend/sprawl_utils/fuzzy_distance.pyx"],
        include_dirs=[np.get_include()] + sprawl_utils_path,
        language="c++",
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
    ),
    Extension(
        "PartSegCore_compiled_backend.color_image_cython",
        ["src/PartSegCore_compiled_backend/color_image_cython.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
        language="c++",
    ),
    Extension(
        "PartSegCore_compiled_backend.multiscale_opening.mso_bind",
        sources=["src/PartSegCore_compiled_backend/multiscale_opening/mso_bind.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
        language="c++",
    ),
]

setup(
    ext_modules=extensions,
    include_package_data=True,
    use_scm_version=True,
)
