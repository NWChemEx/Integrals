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
#include "integrals/uncertain_types.hpp"

using namespace integrals::testing;

using udouble            = integrals::type::uncertain_double;
constexpr bool has_sigma = integrals::type::has_sigma();

template<typename InputType>
auto unwrap_mean(InputType&& uq) {
    if constexpr(has_sigma) {
        return uq.mean();
    } else {
        return uq;
    }
}

template<typename InputType>
auto unwrap_sd(InputType&& uq) {
    if constexpr(has_sigma) {
        return uq.sd();
    } else {
        return 0.0;
    }
}

TEST_CASE("OperatorBase") {
    using aos_t         = simde::type::aos;
    using aos_squared_t = simde::type::aos_squared;
    using op_t          = simde::type::v_ee_type;
    using op_base_t     = simde::type::op_base_type;

    auto mm = initialize_integrals();
    REQUIRE(mm.count("Evaluate 2-Index BraKet"));
    REQUIRE(mm.count("Evaluate 3-Index BraKet"));
    REQUIRE(mm.count("Evaluate 4-Index BraKet"));

    // Get basis set
    auto mol  = water_molecule();
    auto aobs = water_sto3g_basis_set();

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
        auto& mod = mm.at("Evaluate 2-Index BraKet");

        SECTION("No UQ") {
            auto T = mod.run_as<test_pt>(braket);

            // Check output
            REQUIRE(trace<2>(T.buffer()) ==
                    Catch::Approx(124.7011973877891364).margin(1.0e-16));
            REQUIRE(norm<2>(T.buffer()) ==
                    Catch::Approx(90.2562579028763707).margin(1.0e-16));
        }

        SECTION("With UQ") {
            if constexpr(has_sigma) {
                mod.change_input("With UQ?", true);
                auto T = mod.run_as<test_pt>(braket);

                // Check output
                REQUIRE(unwrap_mean(trace<2, udouble>(T.buffer())) ==
                        Catch::Approx(124.7011973877891364).margin(1.0e-16));
                REQUIRE(unwrap_sd(trace<2, udouble>(T.buffer())) ==
                        Catch::Approx(7e-16).margin(1.0e-16));
                REQUIRE(unwrap_mean(norm<2, udouble>(T.buffer())) ==
                        Catch::Approx(90.2562579028763707).margin(1.0e-16));
                REQUIRE(unwrap_sd(norm<2, udouble>(T.buffer())) ==
                        Catch::Approx(3e-16).margin(1.0e-16));
            }
        }
    }
    SECTION("3-Index") {
        using braket_t = simde::type::braket<aos_t, op_base_t, aos_squared_t>;
        using test_pt  = simde::EvaluateBraKet<braket_t>;

        // Make BraKet Input
        braket_t braket(aos, op_base, aos_squared);

        auto& mod = mm.at("Evaluate 3-Index BraKet");

        SECTION("No UQ") {
            // Call module
            auto T = mod.run_as<test_pt>(braket);

            // Check output
            REQUIRE(trace<3>(T.buffer()) ==
                    Catch::Approx(16.8245948391706577).margin(1.0e-16));
            REQUIRE(norm<3>(T.buffer()) ==
                    Catch::Approx(20.6560572032543597).margin(1.0e-16));
        }

        SECTION("With UQ") {
            if constexpr(has_sigma) {
                mod.change_input("With UQ?", true);
                // Call module
                auto T = mod.run_as<test_pt>(braket);

                // Check output
                auto& t = T.buffer();
                REQUIRE(unwrap_mean(trace<3, udouble>(t)) ==
                        Catch::Approx(16.8245948391706577).margin(1.0e-16));
                REQUIRE(unwrap_sd(trace<3, udouble>(t)) ==
                        Catch::Approx(7e-16).margin(1.0e-16));
                REQUIRE(unwrap_mean(norm<3, udouble>(t)) ==
                        Catch::Approx(20.6560572032543597).margin(1.0e-16));
                REQUIRE(unwrap_sd(norm<3, udouble>(t)) ==
                        Catch::Approx(7e-16).margin(1.0e-16));
            }
        }
    }

    SECTION("4-Index") {
        using braket_t =
          simde::type::braket<aos_squared_t, op_base_t, aos_squared_t>;
        using test_pt = simde::EvaluateBraKet<braket_t>;

        // Make BraKet Input
        braket_t braket(aos_squared, op_base, aos_squared);

        auto& mod = mm.at("Evaluate 4-Index BraKet");

        SECTION("No UQ") {
            // Call module
            auto T = mod.run_as<test_pt>(braket);

            // Check output
            auto& t = T.buffer();
            REQUIRE(trace<4>(t) ==
                    Catch::Approx(9.7919608941952063).margin(1.0e-16));
            REQUIRE(norm<4>(t) ==
                    Catch::Approx(7.7796143419802553).margin(1.0e-16));
        }

        SECTION("With UQ") {
            if constexpr(has_sigma) {
                // Call module
                mod.change_input("With UQ?", true);
                auto T = mod.run_as<test_pt>(braket);

                // Check output
                auto& t = T.buffer();
                REQUIRE(unwrap_mean(trace<4, udouble>(t)) ==
                        Catch::Approx(9.7919608941952063).margin(1.0e-16));
                REQUIRE(unwrap_sd(trace<4, udouble>(t)) ==
                        Catch::Approx(7e-16).margin(1.0e-16));
                REQUIRE(unwrap_mean(norm<4, udouble>(t)) ==
                        Catch::Approx(7.7796143419802553).margin(1.0e-16));
                REQUIRE(unwrap_sd(norm<4, udouble>(t)) ==
                        Catch::Approx(11e-16).margin(1.0e-16));
            }
        }
    }
}
