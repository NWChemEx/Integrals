/*
 * Copyright 2025 NWChemEx-Project
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

TEST_CASE("Density Fitting Integral") {
    using test_pt = simde::ERI3;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("Density Fitting Integral"));

    // Get basis set
    auto mol  = test::h2_molecule();
    auto aobs = test::h2_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);
    simde::type::aos_squared aos_squared(aos, aos);

    // Make Operator
    simde::type::v_ee_type op{};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos_squared);

    // Call module
    auto T = mm.at("Density Fitting Integral").run_as<test_pt>(braket);

    // Check output
    // auto t = test::eigen_buffer<3>(T.buffer());
}