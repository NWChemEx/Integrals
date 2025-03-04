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

TEST_CASE("Coulomb Metric") {
    using test_pt = simde::ERI2;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("Coulomb Metric"));

    // Get basis set
    auto mol  = test::h2_molecule();
    auto aobs = test::h2_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Make Operator
    simde::type::v_ee_type op{};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos);

    // Call module
    auto T = mm.at("Coulomb Metric").run_as<test_pt>(braket);

    auto t = test::eigen_tensor<2>(T.buffer());
    REQUIRE(t(0, 0) == Catch::Approx(0.15824726).margin(1E-6));
    REQUIRE(t(0, 1) == Catch::Approx(0.0).margin(1E-6));
    REQUIRE(t(1, 0) == Catch::Approx(-0.23097095).margin(1E-6));
    REQUIRE(t(1, 1) == Catch::Approx(0.27998174).margin(1E-6));
}
