cmake_minimum_required(VERSION 3.8)
project(integrals VERSION 1.0.0)
find_package(cpp REQUIRED)
cpp_option(BUILD_TESTS ON)

cpp_find_or_build_dependency(
    NAME property_types
    URL github.com/NWChemEx-Project/PropertyTypes
    PRIVATE
    CMAKE_ARGS BUILD_TESTS=OFF
)

cpp_find_dependency(
    NAME libint
    FIND_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Findlibint.cmake
)

add_subdirectory(integrals)

### Testing ###
if(BUILD_TESTS)
    cpp_find_or_build_dependency(NAME Catch2 URL github.com/catchorg/Catch2)
    enable_testing()
    add_subdirectory(tests)
endif()
