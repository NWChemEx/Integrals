/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
