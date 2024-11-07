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

#include "../water_sto3g.hpp"
#include <catch2/catch_test_macros.hpp>
#include <integrals/integrals.hpp>

namespace {
simde::type::tensor correct_value() {
    return simde::type::tensor{
      {
        1.0464370899978459,
        3.4291996305312606,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        2.6052624057150817,
        2.6052624057150817,
      },
      {
        3.4291996305312606,
        26.435225216427671,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        25.3420821293274088,
        25.3420821293274088,
      },
      {
        0.0000000000000000,
        0.0000000000000000,
        5.7847978365504318,
        0.0000000000000000,
        0.0000000000000000,
        3.2924421173969143,
        3.2924421173969143,
      },
      {
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        5.7847978365504300,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
      },
      {
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        5.7847978365504318,
        4.2141100538676941,
        -4.2141100538676941,
      },
      {
        2.6052624057150817,
        25.3420821293274088,
        3.2924421173969143,
        0.0000000000000000,
        4.2141100538676941,
        39.9325707858561643,
        26.6712894368540034,
      },
      {
        2.6052624057150817,
        25.3420821293274088,
        3.2924421173969143,
        0.0000000000000000,
        -4.2141100538676941,
        26.6712894368540034,
        39.9325707858561643,
      },
    };
}
} // namespace

TEST_CASE("ERI2") {
    using test_pt = simde::ERI2;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("ERI2"));

    // Get basis set
    auto mol  = test::water_molecule();
    auto aobs = test::water_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Make Operator
    simde::type::v_ee_type op{};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos);

    // Call module
    auto T    = mm.at("ERI2").run_as<test_pt>(braket);
    auto corr = correct_value();

    // Check output
    REQUIRE(T == corr);
}
