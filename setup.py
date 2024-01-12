import os
import platform

import numpy as np
from setuptools import Extension, setup

current_dir = os.path.dirname(os.path.abspath(__file__))
package_dir = os.path.join(current_dir, 'src')

cpp_standard = ['-std=c++17', '-g0', '-O2', '-DNDEBUG']  # "-DDEBUG", "-O0", "-ggdb3" ]
sprawl_utils_path = [os.path.join(package_dir, 'PartSegCore_compiled_backend', 'sprawl_utils')]

if platform.system() == 'Darwin':
    cpp_standard += ['-stdlib=libc++', '-mmacosx-version-min=10.9', '-Wno-nullability-completeness']
    omp = ['-Xpreprocesssor', '-fopenmp']
elif platform.system() == 'Linux':
    omp = ['-fopenmp']
else:
    omp = ["/openmp"]



extensions = [
    Extension(
        'PartSegCore_compiled_backend.sprawl_utils.euclidean_cython',
        sources=['src/PartSegCore_compiled_backend/sprawl_utils/euclidean_cython.pyx'],
        include_dirs=[np.get_include()] + sprawl_utils_path,
        language='c++',
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
    ),
    Extension(
        'PartSegCore_compiled_backend.sprawl_utils.path_sprawl_cython',
        sources=['src/PartSegCore_compiled_backend/sprawl_utils/path_sprawl_cython.pyx'],
        include_dirs=[np.get_include()] + sprawl_utils_path,
        language='c++',
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
    ),
    Extension(
        'PartSegCore_compiled_backend.sprawl_utils.sprawl_utils',
        sources=['src/PartSegCore_compiled_backend/sprawl_utils/sprawl_utils.pyx'],
        include_dirs=[np.get_include()] + sprawl_utils_path,
        language='c++',
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
    ),
    Extension(
        'PartSegCore_compiled_backend.sprawl_utils.fuzzy_distance',
        sources=['src/PartSegCore_compiled_backend/sprawl_utils/fuzzy_distance.pyx'],
        include_dirs=[np.get_include()] + sprawl_utils_path,
        language='c++',
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
    ),
    Extension(
        'PartSegCore_compiled_backend.color_image_cython',
        ['src/PartSegCore_compiled_backend/color_image_cython.pyx'],
        include_dirs=[np.get_include()],
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
        language='c++',
    ),
    Extension(
        'PartSegCore_compiled_backend.utils',
        ['src/PartSegCore_compiled_backend/utils.pyx'],
        include_dirs=[np.get_include()],
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
        language='c++',
    ),
    Extension(
        'PartSegCore_compiled_backend.multiscale_opening.mso_bind',
        sources=['src/PartSegCore_compiled_backend/multiscale_opening/mso_bind.pyx'],
        include_dirs=[np.get_include()],
        extra_compile_args=cpp_standard,
        extra_link_args=cpp_standard,
        language='c++',
    ),
    Extension(
        'PartSegCore_compiled_backend._fast_unique',
        sources=['src/PartSegCore_compiled_backend/_fast_unique.pyx'],
        include_dirs=[np.get_include()],
        extra_compile_args=cpp_standard + omp
        if platform.system() == 'Darwin'
        else ['-fopenmp'],
        extra_link_args=cpp_standard + omp,
        language='c++',
    ),
]

setup(ext_modules=extensions, include_package_data=True, use_scm_version=True)
