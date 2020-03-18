#include "test_common_TA.hpp"
#include <integrals/nwx_direct/eri_direct_type.hpp>
#include <integrals/nwx_direct/stg_direct_type.hpp>
#include <integrals/nwx_direct/yukawa_direct_type.hpp>

TEST_CASE("Direct ERI3") {
    test_property_type<property_types::ERI3CDirect<>>(
            {"Bra", "Ket1", "Ket2", "Derivative"},
            {"ERIs"});
}

TEST_CASE("Direct ERI4") {
    test_property_type<property_types::ERI4CDirect<>>(
            {"Bra1", "Bra2", "Ket1", "Ket2", "Derivative"},
            {"ERIs"});
}

TEST_CASE("Direct STG3") {
    test_property_type<property_types::STG3CDirect<>>(
            {"Bra", "Ket1", "Ket2", "Derivative", "STG Exponent"},
            {"STG Integrals"});
}

TEST_CASE("Direct STG4") {
    test_property_type<property_types::STG4CDirect<>>(
            {"Bra1", "Bra2", "Ket1", "Ket2", "Derivative", "STG Exponent"},
            {"STG Integrals"});
}

TEST_CASE("Direct Yukawa3") {
    test_property_type<property_types::Yukawa3CDirect<>>(
            {"Bra", "Ket1", "Ket2", "Derivative", "STG Exponent"},
            {"Yukawa Integrals"});
}

TEST_CASE("Direct Yukawa4") {
    test_property_type<property_types::Yukawa4CDirect<>>(
            {"Bra1", "Bra2", "Ket1", "Ket2", "Derivative", "STG Exponent"},
            {"Yukawa Integrals"});
}
