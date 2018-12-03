cmake_minimum_required(VERSION 3.8)
project(integrals VERSION 1.0.0)
set(CPP_NAMESPACE Integrals)
find_package(cpp REQUIRED)
cpp_option(BUILD_TESTS ON)
cpp_find_or_build_dependency(NAME Catch2 URL github.com/catchorg/Catch2)
cpp_find_or_build_dependency(
    NAME libint
    URL github.com/evaleev/libint
    VERSION v2.4.2
    FIND_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindLibInt.cmake
)
cpp_find_or_build_dependency(
    NAME libchemist
    URL github.com/NWChemEx-Project/LibChemist
    PRIVATE
    CMAKE_ARGS BUILD_TESTS=OFF
)
add_subdirectory(Integrals)

### Testing ###
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(Integrals_Test)
endif()
