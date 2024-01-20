import os
import platform

import numpy as np
from setuptools import Extension, setup

current_dir = os.path.dirname(os.path.abspath(__file__))
package_dir = os.path.join(current_dir, 'src')

cpp_standard = ['-std=c++11', '-g0', '-O2', '-DNDEBUG']  # "-DDEBUG", "-O0", "-ggdb3" ]
sprawl_utils_path = [os.path.join(package_dir, 'PartSegCore_compiled_backend', 'sprawl_utils')]
extra_link_args = []

if platform.system() == 'Darwin':
    cpp_standard += ['-stdlib=libc++', '-Wno-nullability-completeness']
    omp = ['-Xpreprocessor', '-fopenmp', '-lomp']
    if omp_prefix := os.environ.get('OMP'):
        cpp_standard += ['-I' + os.path.join(omp_prefix, 'include')]
        extra_link_args += ['-L' + os.path.join(omp_prefix, 'lib')]
elif platform.system() == 'Linux':
    omp = ['-fopenmp']
else:
    omp = ['/openmp']


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
        'PartSegCore_compiled_backend.calc_bounds',
        ['src/PartSegCore_compiled_backend/calc_bounds.pyx'],
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
        extra_compile_args=cpp_standard + omp,
        extra_link_args=cpp_standard + omp + extra_link_args,
        language='c++',
    ),
    Extension(
        'PartSegCore_compiled_backend._napari_mapping',
        sources=['src/PartSegCore_compiled_backend/_napari_mapping.pyx'],
        include_dirs=[np.get_include()],
        extra_compile_args=cpp_standard + omp,
        extra_link_args=cpp_standard + omp + extra_link_args,
        language='c++',
    ),
]

setup(ext_modules=extensions, include_package_data=True, use_scm_version=True)
