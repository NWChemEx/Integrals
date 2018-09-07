set(
    INTEGRALS_TESTS test_main.cpp
                    TestCommon.hpp
                    TestDF_libint.cpp
                    TestERI_libint.cpp
                    TestKinetic_libint.cpp
                    TestMetric_libint.cpp
                    TestNuclearElectron_libint.cpp
                    TestOverlap_libint.cpp
                    TestEDipole_libint.cpp
)

add_catch_cxx_tests(
        NAME test_integrals
        SOURCES ${INTEGRALS_TESTS}
        DEPENDS Integrals
)
