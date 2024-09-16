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

# Downloads common CMake modules used throughout NWChemEx
include(cmake/get_nwx_cmake.cmake)


#Sets the version to whatever git thinks it is
include(get_version_from_git)
get_version_from_git(integrals_version "${CMAKE_CURRENT_LIST_DIR}")
project(integrals integrals_version "${VERSION}" LANGUAGES CXX)

include(nwx_versions)
include(get_cmaize)
include(nwx_cxx_api_docs)

### Paths ###
set(INTEGRALS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INTEGRALS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(INTEGRALS_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")

nwx_cxx_api_docs("${INTEGRALS_SOURCE_DIR}" "${INTEGRALS_INCLUDE_DIR}")

### Options ###
cmaize_option_list(
    BUILD_TESTING OFF "Should we build the tests?"
    BUILD_PYBIND11_PYBINDINGS ON "Build pybind11 python3 bindings?"
)

# Can't build from github due to Libint setup.
# Will find prebuilt. Probably for the best.
cmaize_find_or_build_dependency(
    Libint2
    URL github.com/evaleev/libint
    BUILD_TARGET int2
    FIND_TARGET Libint2::int2
)

cmaize_find_or_build_dependency(
    simde
    URL github.com/NWChemEx/SimDE
    BUILD_TARGET simde
    FIND_TARGET nwx::simde
    CMAKE_ARGS BUILD_TESTING=OFF
               BUILD_PYBIND11_PYBINDINGS=${BUILD_PYBIND11_PYBINDINGS}
)

cmaize_add_library(
    ${PROJECT_NAME}
    SOURCE_DIR "${INTEGRALS_SOURCE_DIR}/${PROJECT_NAME}"
    INCLUDE_DIRS "${INTEGRALS_INCLUDE_DIR}/${PROJECT_NAME}"
    DEPENDS Libint2 simde
)

include(nwx_pybind11)
nwx_add_pybind11_module(
    ${PROJECT_NAME}
    SOURCE_DIR "${INTEGRALS_SOURCE_DIR}/python"
    DEPENDS "${PROJECT_NAME}"
)

if("${BUILD_TESTING}")
    set(CXX_TEST_DIR ${INTEGRALS_TESTS_DIR}/tests/cxx/)
    set(PYTHON_TEST_DIR ${INTEGRALS_TESTS_DIR}/tests/python)

    cmaize_find_or_build_dependency(
        Catch2
        URL github.com/catchorg/Catch2
        BUILD_TARGET Catch2
        FIND_TARGET Catch2::Catch2
        VERSION v3.6.0
    )

    cmaize_add_tests(
        test_unit_integrals
        SOURCE_DIR "${CXX_TEST_DIR}/unit"
        INCLUDE_DIRS "${INTEGRALS_SOURCE_DIR}/${PROJECT_NAME}"
        DEPENDS Catch2 ${PROJECT_NAME}
    )
    
    nwx_pybind11_tests(
        py_integrals "${PYTHON_TEST_DIR}/unit_tests/test_integrals.py"
        SUBMODULES pluginplay chemist parallelzone simde
    )
endif()

cmaize_add_package(${PROJECT_NAME} NAMESPACE nwx::)
