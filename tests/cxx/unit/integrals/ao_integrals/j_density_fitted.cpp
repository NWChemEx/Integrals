/*
 * Copyright 2024 NWChemEx-Project
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

TEST_CASE("Density Fitted J builder") {
    using pt = simde::aos_j_e_aos;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("Density Fitted J builder"));

    // Get basis set
    auto mol  = test::h2_molecule();
    auto aobs = test::h2_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Make Operator
    simde::type::j_e_type op(simde::type::electron{}, test::h2_density());

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos);

    // Call module
    mm.change_input("Density Fitted J builder", "Auxiliary Basis Set", aos);
    const auto& T = mm.at("Density Fitted J builder").run_as<pt>(braket);

    auto t = test::eigen_tensor<2>(T.buffer());
    REQUIRE(t(0, 0) == Catch::Approx(0.50515668).margin(1E-6));
    REQUIRE(t(0, 1) == Catch::Approx(0.22344536).margin(1E-6));
    REQUIRE(t(1, 0) == Catch::Approx(0.22344536).margin(1E-6));
    REQUIRE(t(1, 1) == Catch::Approx(0.50515668).margin(1E-6));
}