set(
    INTEGRALS_TESTS test_main.cpp
                    TestCommon.hpp
                    TestOverlap_libint.cpp
)

add_catch_cxx_tests(
        NAME test_integrals
        SOURCES ${INTEGRALS_TESTS}
        DEPENDS Integrals
)
