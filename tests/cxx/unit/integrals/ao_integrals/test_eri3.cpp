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

#include "test_ao_integrals.hpp"

TEST_CASE("ERI3") {
    using test_pt = simde::ERI3;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("ERI3"));

    // Get basis set
    auto mol  = test::water_molecule();
    auto aobs = test::water_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);
    simde::type::aos_squared aos_squared(aos, aos);

    // Make Operator
    simde::type::v_ee_type op{};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos_squared);

    // Call module
    auto T = mm.at("ERI3").run_as<test_pt>(braket);

    // Check output
    auto t = test::eigen_buffer<3>(T.buffer());
    REQUIRE(test::trace(t) ==
            Catch::Approx(16.8245948391706577).margin(1.0e-16));
    REQUIRE(test::norm(t) ==
            Catch::Approx(20.6560572032543597).margin(1.0e-16));
}
