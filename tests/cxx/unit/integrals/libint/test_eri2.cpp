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

using namespace integrals;

TEST_CASE("ERI2") {
    using test_pt = simde::ERI2;

    auto mm = integrals::testing::initialize_integrals();
    REQUIRE(mm.count("ERI2"));

    // Get basis set
    auto mol  = testing::water_molecule();
    auto aobs = testing::water_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Make Operator
    simde::type::v_ee_type op{};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos);

    // Call module
    auto T = mm.at("ERI2").run_as<test_pt>(braket);

    // Check output
    auto& t = T.buffer();
    REQUIRE(testing::trace<2>(t) ==
            Catch::Approx(124.7011973877891364).margin(1.0e-16));
    REQUIRE(testing::norm<2>(t) ==
            Catch::Approx(90.2562579028763707).margin(1.0e-16));
}
