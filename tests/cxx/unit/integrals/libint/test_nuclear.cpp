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

#include "../testing/testing.hpp"

using namespace integrals::testing;

TEST_CASE("Nuclear") {
    using test_pt = simde::aos_v_en_aos;

    auto mm = initialize_integrals();
    REQUIRE(mm.count("Nuclear"));

    // Get basis set
    auto mol  = water_molecule();
    auto aobs = water_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Make Operator
    simde::type::v_en_type op{chemist::Electron{}, mol.nuclei().as_nuclei()};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos);

    // Call module
    auto T = mm.at("Nuclear").run_as<test_pt>(braket);

    // Check output
    auto& t = T.buffer();
    REQUIRE(trace<2>(t) ==
            Catch::Approx(-111.9975421879705664).margin(1.0e-16));
    REQUIRE(norm<2>(t) == Catch::Approx(66.4857539908047528).margin(1.0e-16));
}
