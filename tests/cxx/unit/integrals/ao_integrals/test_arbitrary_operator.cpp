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

TEST_CASE("OperatorBase") {
    using aos_t         = simde::type::aos;
    using aos_squared_t = simde::type::aos_squared;
    using op_t          = simde::type::v_ee_type;
    using op_base_t     = simde::type::op_base_type;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("Evaluate 2-Index BraKet"));
    REQUIRE(mm.count("Evaluate 3-Index BraKet"));
    REQUIRE(mm.count("Evaluate 4-Index BraKet"));

    // Get basis set
    auto mol  = test::water_molecule();
    auto aobs = test::water_sto3g_basis_set();

    // Make AOS object
    aos_t aos(aobs);
    aos_squared_t aos_squared(aos, aos);

    // Make Operator
    op_t op{};
    op_base_t& op_base = op;

    SECTION("2-Index") {
        using braket_t = simde::type::braket<aos_t, op_base_t, aos_t>;
        using test_pt  = simde::EvaluateBraKet<braket_t>;

        // Make BraKet Input
        braket_t braket(aos, op_base, aos);

        // Call module
        auto T = mm.at("Evaluate 2-Index BraKet").run_as<test_pt>(braket);

        // Check output
        auto t = test::eigen_buffer<2>(T.buffer());
        REQUIRE(test::trace(t) ==
                Catch::Approx(124.7011973877891364).margin(1.0e-16));
        REQUIRE(test::norm(t) ==
                Catch::Approx(90.2562579028763707).margin(1.0e-16));
    }

    SECTION("3-Index") {
        using braket_t = simde::type::braket<aos_t, op_base_t, aos_squared_t>;
        using test_pt  = simde::EvaluateBraKet<braket_t>;

        // Make BraKet Input
        braket_t braket(aos, op_base, aos_squared);

        // Call module
        auto T = mm.at("Evaluate 3-Index BraKet").run_as<test_pt>(braket);

        // Check output
        auto t = test::eigen_buffer<3>(T.buffer());
        REQUIRE(test::trace(t) ==
                Catch::Approx(16.8245948391706577).margin(1.0e-16));
        REQUIRE(test::norm(t) ==
                Catch::Approx(20.6560572032543597).margin(1.0e-16));
    }

    SECTION("4-Index") {
        using braket_t =
          simde::type::braket<aos_squared_t, op_base_t, aos_squared_t>;
        using test_pt = simde::EvaluateBraKet<braket_t>;

        // Make BraKet Input
        braket_t braket(aos_squared, op_base, aos_squared);

        // Call module
        auto T = mm.at("Evaluate 4-Index BraKet").run_as<test_pt>(braket);

        // Check output
        auto t = test::eigen_buffer<4>(T.buffer());
        REQUIRE(test::trace(t) ==
                Catch::Approx(9.7919608941952063).margin(1.0e-16));
        REQUIRE(test::norm(t) ==
                Catch::Approx(7.7796143419802553).margin(1.0e-16));
    }
}
