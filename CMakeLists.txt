cmake_minimum_required(VERSION 3.21)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-DDEBUG -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG")


# Define compiler directive
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

if(DEFINED ENV{CUSTOM_MODULE_PATH})
    list(APPEND CMAKE_MODULE_PATH $ENV{CUSTOM_MODULE_PATH})
    message("Custom module path: ${CMAKE_MODULE_PATH}")
endif()


find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(Cython MODULE REQUIRED VERSION 3.0)
include(UseCython)

if(DEFINED ENV{OMP})
    list(APPEND CMAKE_PREFIX_PATH $ENV{OMP})
endif()

# find OpenMP
find_package(OpenMP REQUIRED)


execute_process(
        COMMAND "${Python_EXECUTABLE}"
        -c "import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)


cython_transpile(src/PartSegCore_compiled_backend/sprawl_utils/euclidean_cython.pyx LANGUAGE CXX OUTPUT_VARIABLE euclidean_cython_cxx)
python_add_library(euclidean_cython MODULE
        "${euclidean_cython_cxx}"
        src/PartSegCore_compiled_backend/sprawl_utils/my_queue.h
        src/PartSegCore_compiled_backend/sprawl_utils/global_consts.h
        WITH_SOABI)
target_include_directories(euclidean_cython PRIVATE "${NUMPY_INCLUDE_DIR}")
target_include_directories(euclidean_cython PRIVATE src/PartSegCore_compiled_backend/sprawl_utils/)

cython_transpile(src/PartSegCore_compiled_backend/sprawl_utils/path_sprawl_cython.pyx LANGUAGE CXX OUTPUT_VARIABLE path_sprawl_cython_cxx)
python_add_library(path_sprawl_cython MODULE
        "${path_sprawl_cython_cxx}"
        src/PartSegCore_compiled_backend/sprawl_utils/my_queue.h
        src/PartSegCore_compiled_backend/sprawl_utils/global_consts.h
        WITH_SOABI)
target_include_directories(path_sprawl_cython PRIVATE "${NUMPY_INCLUDE_DIR}")
target_include_directories(path_sprawl_cython PRIVATE src/PartSegCore_compiled_backend/sprawl_utils/)

cython_transpile(src/PartSegCore_compiled_backend/sprawl_utils/sprawl_utils.pyx LANGUAGE CXX OUTPUT_VARIABLE sprawl_utils_cython_cxx)
python_add_library(sprawl_utils MODULE
        "${sprawl_utils_cython_cxx}"
        src/PartSegCore_compiled_backend/sprawl_utils/my_queue.h
        src/PartSegCore_compiled_backend/sprawl_utils/global_consts.h
        WITH_SOABI)
target_include_directories(sprawl_utils PRIVATE "${NUMPY_INCLUDE_DIR}")
target_include_directories(sprawl_utils PRIVATE src/PartSegCore_compiled_backend/sprawl_utils/)


cython_transpile(src/PartSegCore_compiled_backend/sprawl_utils/fuzzy_distance.pyx LANGUAGE CXX OUTPUT_VARIABLE fuzzy_distance_cython_cxx)
python_add_library(fuzzy_distance MODULE
        "${fuzzy_distance_cython_cxx}"
        src/PartSegCore_compiled_backend/sprawl_utils/my_queue.h
        src/PartSegCore_compiled_backend/sprawl_utils/global_consts.h
        WITH_SOABI)
target_include_directories(fuzzy_distance PRIVATE "${NUMPY_INCLUDE_DIR}")
target_include_directories(fuzzy_distance PRIVATE src/PartSegCore_compiled_backend/sprawl_utils/)


cython_transpile(src/PartSegCore_compiled_backend/color_image_cython.pyx LANGUAGE CXX OUTPUT_VARIABLE color_image_cython_cxx)
python_add_library(color_image_cython MODULE "${color_image_cython_cxx}" WITH_SOABI)
target_include_directories(color_image_cython PRIVATE "${NUMPY_INCLUDE_DIR}")


cython_transpile(src/PartSegCore_compiled_backend/calc_bounds.pyx LANGUAGE CXX OUTPUT_VARIABLE calc_bounds_cython_cxx)
python_add_library(calc_bounds MODULE "${calc_bounds_cython_cxx}" WITH_SOABI)
target_include_directories(calc_bounds PRIVATE "${NUMPY_INCLUDE_DIR}")


cython_transpile(src/PartSegCore_compiled_backend/multiscale_opening/mso_bind.pyx LANGUAGE CXX OUTPUT_VARIABLE mso_bind_cython_cxx)
python_add_library(
    mso_bind MODULE
    "${mso_bind_cython_cxx}"
    # src/PartSegCore_compiled_backend/multiscale_opening/mso.cpp
    src/PartSegCore_compiled_backend/multiscale_opening/mso.h
    src/PartSegCore_compiled_backend/multiscale_opening/my_queue.h
    WITH_SOABI)
target_include_directories(mso_bind PRIVATE src/PartSegCore_compiled_backend/multiscale_opening/)
target_include_directories(mso_bind PRIVATE "${NUMPY_INCLUDE_DIR}")


cython_transpile(src/PartSegCore_compiled_backend/_fast_unique.pyx LANGUAGE CXX OUTPUT_VARIABLE fast_unique_cython_cxx)
python_add_library(_fast_unique MODULE "${fast_unique_cython_cxx}" WITH_SOABI)
target_include_directories(_fast_unique PRIVATE "${NUMPY_INCLUDE_DIR}")
target_link_libraries(_fast_unique PRIVATE OpenMP::OpenMP_CXX)


cython_transpile(src/PartSegCore_compiled_backend/_napari_mapping.pyx LANGUAGE CXX OUTPUT_VARIABLE napari_mapping_cython_cxx)
python_add_library(_napari_mapping MODULE "${napari_mapping_cython_cxx}" WITH_SOABI)
target_include_directories(_napari_mapping PRIVATE "${NUMPY_INCLUDE_DIR}")
target_link_libraries(_napari_mapping PRIVATE OpenMP::OpenMP_CXX)

cython_transpile(src/PartSegCore_compiled_backend/triangulate.pyx LANGUAGE CXX OUTPUT_VARIABLE triangulate_cython_cxx)
python_add_library(
        triangulate MODULE
        "${triangulate_cython_cxx}"
        src/PartSegCore_compiled_backend/triangulation/triangulate.hpp
        src/PartSegCore_compiled_backend/triangulation/intersection.hpp
        src/PartSegCore_compiled_backend/triangulation/point.hpp
        src/PartSegCore_compiled_backend/triangulation/debug_util.hpp
        WITH_SOABI)
target_include_directories(triangulate PRIVATE "${NUMPY_INCLUDE_DIR}")
target_include_directories(triangulate PRIVATE src/PartSegCore_compiled_backend/)


install(TARGETS euclidean_cython DESTINATION PartSegCore_compiled_backend/sprawl_utils)
install(TARGETS path_sprawl_cython DESTINATION PartSegCore_compiled_backend/sprawl_utils)
install(TARGETS sprawl_utils DESTINATION PartSegCore_compiled_backend/sprawl_utils)
install(TARGETS fuzzy_distance DESTINATION PartSegCore_compiled_backend/sprawl_utils)
install(TARGETS color_image_cython DESTINATION PartSegCore_compiled_backend)
install(TARGETS calc_bounds DESTINATION PartSegCore_compiled_backend)
install(TARGETS mso_bind DESTINATION PartSegCore_compiled_backend/multiscale_opening)
install(TARGETS _fast_unique DESTINATION PartSegCore_compiled_backend)
install(TARGETS _napari_mapping DESTINATION PartSegCore_compiled_backend)
install(TARGETS triangulate DESTINATION PartSegCore_compiled_backend)
