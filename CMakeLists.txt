cmake_minimum_required(VERSION 3.14)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" VERSION)
project(integrals VERSION "${VERSION}" LANGUAGES CXX)

set(
        CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
        CACHE STRING "" FORCE
)

include(get_cpp)
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("${CMAKE_CURRENT_SOURCE_DIR}/integrals")

# ### Options ###
# option(BUILD_TESTING "Should we build the tests?" OFF)
# option(BUILD_PYBINDINGS "Should we build Python3 bindings?" ON)
#
#
# # Can't build from github due to Libint setup.
# # Will find prebuilt. Probably for the best.
# cpp_find_or_build_dependency(
#         Libint2
#         URL github.com/evaleev/libint
#         BUILD_TARGET int2
#         FIND_TARGET Libint2::int2
# )
#
# cpp_find_or_build_dependency(
#         simde
#         URL github.com/NWChemEx-Project/SimDE
#         PRIVATE TRUE
#         BUILD_TARGET simde
#         FIND_TARGET nwx::simde
#         CMAKE_ARGS BUILD_TESTING=OFF
# )
#
# set(project_src_dir "${CMAKE_CURRENT_LIST_DIR}/src/${PROJECT_NAME}")
# set(project_inc_dir "${CMAKE_CURRENT_LIST_DIR}/include/${PROJECT_NAME}")
#
# cpp_add_library(
#         ${PROJECT_NAME}
#         SOURCE_DIR "${project_src_dir}"
#         INCLUDE_DIR "${project_inc_dir}"
#         DEPENDS Libint2::int2 nwx::simde
# )
#
# include(nwx_python_mods)
# cppyy_make_python_package(PACKAGE integrals NAMESPACES integrals DEPPACKAGES simde )
#
# if("${BUILD_TESTING}")
#     cpp_find_or_build_dependency(
#             Catch2
#             URL github.com/catchorg/Catch2
#             BUILD_TARGET Catch2
#             VERSION v2.x
#             FIND_TARGET Catch2::Catch2
#     )
#     cpp_find_or_build_dependency(
#         mokup
#         URL github.com/NWChemEx-Project/Mokup
#         PRIVATE TRUE
#         BUILD_TARGET mokup
#         FIND_TARGET nwx::mokup
#         CMAKE_ARGS BUILD_TESTING=OFF
#     )
#     cpp_add_tests(
#         test_integrals
#         SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/tests"
#         INCLUDE_DIR "${project_src_dir}"
#         DEPENDS Catch2::Catch2 nwx::mokup ${PROJECT_NAME}
#     )
# endif()
