set(
    INTEGRALS_TESTS test_main.cpp
#                    TestDF_libint.cpp
                    TestERI_libint.cpp
#                    TestKinetic_libint.cpp
#                    TestMetric_libint.cpp
#                    TestNuclearElectron_libint.cpp
                    TestOverlap_libint.cpp
#                    TestEDipole_libint.cpp
#                    TestEQuadrupole_libint.cpp
#                    TestEOctopole_libint.cpp
)

cpp_add_executable(
    test_integrals
    SOURCES ${INTEGRALS_TESTS}
    DEPENDS integrals Catch2::Catch2
)

add_test(NAME test_integrals COMMAND test_integrals)
