# Copyright 2022 NWChemEx-Project
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
set(VERSION 1.0.0) #TODO: get from git
project(integrals VERSION "${VERSION}" LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
    nwx_cmake
    GIT_REPOSITORY https://github.com/NWChemEx-Project/NWXCMake
)
FetchContent_MakeAvailable(nwx_cmake)
list(APPEND CMAKE_MODULE_PATH "${nwx_cmake_SOURCE_DIR}/cmake")

set(
    CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
    CACHE STRING "" FORCE
)

include(get_cmaize)
include(nwx_versions)
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("${CMAKE_CURRENT_SOURCE_DIR}/integrals")

### Options ###
option(BUILD_TESTING "Should we build the tests?" OFF)
option(BUILD_PYBIND11_PYBINDINGS "Use pybind11 to build Python3 bindings?" OFF)
option(BUILD_CPPYY_PYBINDINGS "Use Cppyy to build Python3 bindings?" OFF)


# Can't build from github due to Libint setup.
# Will find prebuilt. Probably for the best.
cmaize_find_or_build_dependency(
    Libint2
    URL github.com/evaleev/libint
    VERSION ${NWX_LIBINT2_VERSION}
    BUILD_TARGET int2
    FIND_TARGET Libint2::int2
)

cmaize_find_or_build_dependency(
    simde
    URL github.com/NWChemEx-Project/SimDE
    PRIVATE TRUE
    VERSION ${NWX_SIMDE_VERSION}
    BUILD_TARGET simde
    FIND_TARGET nwx::simde
    CMAKE_ARGS BUILD_TESTING=OFF
)

set(project_src_dir "${CMAKE_CURRENT_LIST_DIR}/src/${PROJECT_NAME}")
set(project_inc_dir "${CMAKE_CURRENT_LIST_DIR}/include/${PROJECT_NAME}")

cmaize_add_library(
    ${PROJECT_NAME}
    SOURCE_DIR "${project_src_dir}"
    INCLUDE_DIRS "${project_inc_dir}"
    DEPENDS Libint2 simde
)

include(nwx_pybind11)
nwx_add_pybind11_module(
    ${PROJECT_NAME}
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/python"
    DEPENDS "${PROJECT_NAME}"
)

include(nwx_python_mods)
cppyy_make_python_package(PACKAGE integrals NAMESPACES integrals DEPPACKAGES simde )

if("${BUILD_TESTING}")
    set(CXX_TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/tests/cxx)
    set(PYTHON_TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/tests/python)

    cmaize_find_or_build_dependency(
        Catch2
        URL github.com/catchorg/Catch2
        BUILD_TARGET Catch2
        FIND_TARGET Catch2::Catch2
        VERSION ${NWX_CATCH2_VERSION}
    )
    cmaize_find_or_build_dependency(
        mokup
        URL github.com/NWChemEx-Project/Mokup
        PRIVATE TRUE
        VERSION ${NWX_MOKUP_VERSION}
        BUILD_TARGET mokup
        FIND_TARGET nwx::mokup
        CMAKE_ARGS BUILD_TESTING=OFF
    )
    cmaize_add_tests(
        test_unit_integrals
        SOURCE_DIR "${CXX_TEST_DIR}/unit"
        INCLUDE_DIRS "${project_src_dir}"
        DEPENDS Catch2 mokup ${PROJECT_NAME}
    )
    cmaize_add_tests(
        test_integration_integrals
        SOURCE_DIR "${CXX_TEST_DIR}/integration"
        INCLUDE_DIRS "${project_src_dir}"
        DEPENDS Catch2 mokup ${PROJECT_NAME}
    )
    nwx_pybind11_tests(
        py_integrals "${PYTHON_TEST_DIR}/unit_tests/test_integrals.py"
        SUBMODULES pluginplay chemist parallelzone simde
    )
endif()

cmaize_add_package(${PROJECT_NAME} NAMESPACE nwx::)
