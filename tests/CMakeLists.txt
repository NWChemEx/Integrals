set(
    INTEGRALS_TESTS test_main.cpp
                    edipole.cpp
                    eoctopole.cpp
                    equadrupole.cpp
                    eri3.cpp
                    eri4.cpp
                    kinetic.cpp
                    metric.cpp
                    nuclear_electron.cpp
                    overlap.cpp
)

cpp_add_executable(
    test_integrals
    SOURCES ${INTEGRALS_TESTS}
    DEPENDS integrals Catch2::Catch2
)

add_test(NAME test_integrals COMMAND test_integrals)
