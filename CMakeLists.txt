cmake_minimum_required(VERSION 3.14)
project(integrals VERSION 1.0.0)

set(
        CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
        CACHE STRING "" FORCE
)

include(get_cpp)

### Options ###
option(BUILD_TESTING "Should we build the tests?" OFF)
option(BUILD_PYBINDINGS "Should we build Python3 bindings?" ON)


# Can't build from github due to Libint setup.
# Will find prebuilt. Probably for the best.
cpp_find_or_build_dependency(
        Libint2
        URL github.com/evaleev/libint
        BUILD_TARGET int2
        FIND_TARGET Libint2::int2
)

cpp_find_or_build_dependency(
        property_types
        URL github.com/NWChemEx-Project/PropertyTypes
        PRIVATE TRUE
        BUILD_TARGET property_types
        FIND_TARGET nwx::property_types
        CMAKE_ARGS BUILD_TESTING=OFF
)

cpp_add_library(
        integrals
        SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/integrals"
        INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/integrals"
        DEPENDS Libint2::int2 nwx::property_types
)

if("${BUILD_TESTING}")
    cpp_find_or_build_dependency(
            Catch2
            URL github.com/catchorg/Catch2
            BUILD_TARGET Catch2
            FIND_TARGET Catch2::Catch2
    )
    cpp_add_tests(
            test_integrals
            SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/tests"
            INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/integrals"
            DEPENDS Catch2::Catch2 integrals
    )
endif()