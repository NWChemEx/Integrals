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
        -61.5805952694322,
        -7.410821856331163,
        -0.0144738837457361,
        0.0000000000000000,
        0.0000000000000000,
        -1.231685572142488,
        -1.231685572142488,
      },
      {
        -7.410821856331163,
        -10.00907114207003,
        -0.1768908347336431,
        0.0000000000000000,
        0.0000000000000000,
        -2.977226853578134,
        -2.977226853578134,
      },
      {
        -0.01447388374573611,
        -0.1768908347336431,
        -9.944043341698766,
        0.0000000000000000,
        0.0000000000000000,
        -1.471793338712961,
        -1.471793338712961,
      },
      {
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        -9.875875995090944,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
      },
      {
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        -9.987549935088563,
        -1.822236913476131,
        1.822236913476131,
      },
      {
        -1.231685572142488,
        -2.977226853578134,
        -1.471793338712961,
        0.0000000000000000,
        -1.822236913476131,
        -5.300203252295022,
        -1.067171080472437,
      },
      {
        -1.231685572142488,
        -2.977226853578134,
        -1.471793338712961,
        0.0000000000000000,
        1.82223691347613,
        -1.067171080472437,
        -5.30020325229502,
      },
    };
}
} // namespace

TEST_CASE("Nuclear") {
    using test_pt = simde::aos_v_en_aos;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("Nuclear"));

    // Get basis set
    auto mol  = test::water_molecule();
    auto aobs = test::water_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Make Operator
    simde::type::v_en_type op{chemist::Electron{}, mol.nuclei().as_nuclei()};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos);

    // Call module
    auto T    = mm.at("Nuclear").run_as<test_pt>(braket);
    auto corr = correct_value();

    // Check output
    REQUIRE(T == corr);
}
