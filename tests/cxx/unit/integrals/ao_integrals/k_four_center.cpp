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

#include "../testing/testing.hpp"

using namespace integrals::testing;

TEST_CASE("Four center K builder") {
    using pt = simde::aos_k_e_aos;

    auto mm = initialize_integrals();
    REQUIRE(mm.count("Four center K builder"));

    // Get basis set
    auto mol  = h2_molecule();
    auto aobs = h2_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    std::vector<double> corr{0.45617623, 0.35130947, 0.35130947, 0.45617623};
    SECTION("No QCUP") {
        // Make Operator
        simde::type::k_e_type op(simde::type::electron{}, h2_density<double>());

        // Make BraKet Input
        chemist::braket::BraKet braket(aos, op, aos);

        // Call module
        const auto& T = mm.at("Four center K builder").run_as<pt>(braket);

        auto t = eigen_tensor<2>(T.buffer());
        REQUIRE(t(0, 0) == Catch::Approx(corr[0]).margin(1E-6));
        REQUIRE(t(0, 1) == Catch::Approx(corr[1]).margin(1E-6));
        REQUIRE(t(1, 0) == Catch::Approx(corr[2]).margin(1E-6));
        REQUIRE(t(1, 1) == Catch::Approx(corr[3]).margin(1E-6));
    }

#ifdef ENABLE_SIGMA
    SECTION("QCUP") {
        mm.change_submod("Four center K builder", "Four-center ERI",
                         "UQ Driver");
        mm.change_input("ERI4", "Threshold", 1e-6);

        using uq_type = tensorwrapper::types::idouble;

        // Make Operator
        simde::type::k_e_type op(simde::type::electron{},
                                 h2_density<uq_type>());

        // Make BraKet Input
        chemist::braket::BraKet braket(aos, op, aos);

        // Call module
        const auto& T = mm.at("Four center K builder").run_as<pt>(braket);

        auto t = eigen_tensor<2, uq_type>(T.buffer());

        REQUIRE(t(0, 0).contains(corr[0]));
        REQUIRE(t(0, 1).contains(corr[1]));
        REQUIRE(t(1, 0).contains(corr[2]));
        REQUIRE(t(1, 1).contains(corr[3]));
    }
#endif
}
