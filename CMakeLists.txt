# Copyright 2026 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.14)

# Downloads common CMake modules used throughout NWChemEx. Included before
# project() so nwx_set_version() is available to compute the version.
include(cmake/get_nwx_cmake.cmake)

# Version from scikit-build-core (wheels/sdists) or the latest git tag.
include(nwx_set_version)
nwx_set_version(integrals_version "${CMAKE_CURRENT_LIST_DIR}")
project(integrals VERSION "${integrals_version}" LANGUAGES CXX)

include(disable_in_source_builds)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

option(BUILD_TESTING           "Whether to build the unit tests"    OFF)
option(BUILD_PYBIND11_BINDINGS "Build Python bindings via pybind11" ON)
option(INTEGRALS_ENABLE_OPENMP "Enable OpenMP for threading?"       OFF)

# Documentation target (no-op unless BUILD_DOCS is on)
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("cxx/include" "cxx/src")

### Dependencies ###
include(get_dependencies)
get_dependencies(simde libint2)

set(integrals_private_depends libint2)
if(INTEGRALS_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    list(APPEND integrals_private_depends OpenMP::OpenMP_CXX)
endif()

### Library ###
# simde is a PUBLIC dep (appears in integrals' public headers); libint2 and
# OpenMP are implementation details (PRIVATE).
include(nwx_library)
nwx_library(${PROJECT_NAME} "cxx/include" "cxx/src"
    simde PRIVATE ${integrals_private_depends})

### Python bindings (no-op unless BUILD_PYBIND11_BINDINGS is on) ###
include(nwx_python_module)
nwx_python_module(${PROJECT_NAME} "cxx/src")

### Tests ###
if(BUILD_TESTING)
    include(catch2_tests_from_dir)
    catch2_tests_from_dir(
        "test_unit_${PROJECT_NAME}"
        "tests/cxx/unit/${PROJECT_NAME}"
        ${PROJECT_NAME} ${NWX_DEP_TARGET_libint2}
        PRIVATE_INCLUDES "cxx/src"
    )

    # Demonstration executable exercised alongside the tests.
    add_executable(primitive_error_models
        examples/primitive_error_models/main.cpp
    )
    target_link_libraries(primitive_error_models PRIVATE ${PROJECT_NAME})

    include(nwx_python_test)
    nwx_python_test(
        py_unit_test_${PROJECT_NAME}
        "${CMAKE_CURRENT_LIST_DIR}/tests/python/unit_tests/test_integrals.py"
    )
endif()
