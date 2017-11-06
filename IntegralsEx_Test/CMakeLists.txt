cmake_minimum_required(VERSION 3.1)
project(IntegralsEx-Test VERSION 0.0.0 LANGUAGES CXX)
find_package(IntegralsEx REQUIRED)
include(TargetMacros)

add_cxx_unit_test(nwx_libint/TestOverlap IntegralsEx)
add_cxx_unit_test(nwx_libint/TestKinetic IntegralsEx)
add_cxx_unit_test(nwx_libint/TestMetric IntegralsEx)
add_cxx_unit_test(nwx_libint/TestNuclearElectron IntegralsEx)
add_cxx_unit_test(nwx_libint/TestDF IntegralsEx)
add_cxx_unit_test(nwx_libint/TestERI IntegralsEx)

