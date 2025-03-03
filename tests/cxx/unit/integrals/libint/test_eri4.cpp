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

#include "../testing.hpp"

TEST_CASE("ERI4") {
    using test_pt = simde::ERI4;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("ERI4"));

    // Get basis set
    auto mol  = test::water_molecule();
    auto aobs = test::water_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);
    simde::type::aos_squared aos_squared(aos, aos);

    // Make Operator
    simde::type::v_ee_type op{};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos_squared, op, aos_squared);

    // Call module
    auto T = mm.at("ERI4").run_as<test_pt>(braket);

    // Check output
    auto& t = T.buffer();
    REQUIRE(test::trace<4>(t) ==
            Catch::Approx(9.7919608941952063).margin(1.0e-16));
    REQUIRE(test::norm<4>(t) ==
            Catch::Approx(7.7796143419802553).margin(1.0e-16));
}
