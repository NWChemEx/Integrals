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

    SECTION("No QCUP") {
        // Make Operator
        simde::type::k_e_type op(simde::type::electron{}, h2_density<double>());

        // Make BraKet Input
        chemist::braket::BraKet braket(aos, op, aos);

        // Call module
        const auto& T = mm.at("Four center K builder").run_as<pt>(braket);

        auto t = eigen_tensor<2>(T.buffer());
        REQUIRE(t(0, 0) == Catch::Approx(0.45617623).margin(1E-6));
        REQUIRE(t(0, 1) == Catch::Approx(0.35130947).margin(1E-6));
        REQUIRE(t(1, 0) == Catch::Approx(0.35130947).margin(1E-6));
        REQUIRE(t(1, 1) == Catch::Approx(0.45617623).margin(1E-6));
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

        // Disclaimer: these values are just what was output by the first run
        // they may not actually be correct. FWIW the means are right.
        std::vector<uq_type> corr{
          uq_type{0.456171 - 5.45552e-06, 0.456171 + 5.45552e-06},
          uq_type{0.351305 - 5.45552e-06, 0.351305 + 5.45552e-06},
          uq_type{0.351305 - 5.45552e-06, 0.351305 + 5.45552e-06},
          uq_type{0.456171 - 5.45552e-06, 0.456171 + 5.45552e-06}};

        REQUIRE(t(0, 0).lower() == Catch::Approx(corr[0].lower()).margin(1E-6));
        REQUIRE(t(0, 0).upper() == Catch::Approx(corr[0].upper()).margin(1E-6));
        REQUIRE(t(0, 1).lower() == Catch::Approx(corr[1].lower()).margin(1E-6));
        REQUIRE(t(0, 1).upper() == Catch::Approx(corr[1].upper()).margin(1E-6));
        REQUIRE(t(1, 0).lower() == Catch::Approx(corr[2].lower()).margin(1E-6));
        REQUIRE(t(1, 0).upper() == Catch::Approx(corr[2].upper()).margin(1E-6));
        REQUIRE(t(1, 1).lower() == Catch::Approx(corr[3].lower()).margin(1E-6));
        REQUIRE(t(1, 1).upper() == Catch::Approx(corr[3].upper()).margin(1E-6));
    }
#endif
}
