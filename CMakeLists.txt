cmake_minimum_required(VERSION 3.8)
project(integrals VERSION 1.0.0)
set(CPP_NAMESPACE Integrals)
find_package(cpp REQUIRED)
cpp_option(BUILD_TESTS ON)

cpp_find_or_build_dependency(
    NAME libchemist
    URL github.com/NWChemEx-Project/LibChemist
    BRANCH master
    PRIVATE
    CMAKE_ARGS BUILD_TESTS=OFF
)

cpp_find_dependency(
        NAME libint
        FIND_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Findlibint.cmake
)

add_subdirectory(Integrals)

### Testing ###
if(BUILD_TESTS)
    cpp_find_or_build_dependency(NAME Catch2 URL github.com/catchorg/Catch2)
    enable_testing()
    add_subdirectory(Integrals_Test)
endif()
