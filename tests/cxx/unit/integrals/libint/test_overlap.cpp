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

TEST_CASE("Overlap") {
    using test_pt = simde::aos_s_e_aos;

    pluginplay::ModuleManager mm;
    load_modules(mm);
    REQUIRE(mm.count("Overlap"));

    // Get basis set
    auto mol  = testing::water_molecule();
    auto aobs = testing::water_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Make Operator
    simde::type::s_e_type op{};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos);

    // Call module
    auto S = mm.at("Overlap").run_as<test_pt>(braket);

    // Check output
    auto& t = S.buffer();
    REQUIRE(testing::trace<2>(t) ==
            Catch::Approx(7.00000000000000266).margin(1.0e-16));
    REQUIRE(testing::norm<2>(t) ==
            Catch::Approx(2.87134497074907324).margin(1.0e-16));
}
