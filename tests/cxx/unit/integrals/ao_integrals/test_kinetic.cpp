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

namespace test {
simde::type::tensor correct_value() {
    return simde::type::tensor{
      {
        29.00319994553958,
        -0.1680109393164922,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        -0.008416385185447427,
        -0.008416385185447427,
      },
      {
        -0.1680109393164923,
        0.8081279549303477,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        0.07051733851899882,
        0.07051733851899882,
      },
      {
        0.0000000000000000,
        0.0000000000000000,
        2.528731198194765,
        0.0000000000000000,
        0.0000000000000000,
        0.1149203802569082,
        0.1149203802569082,
      },
      {
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        2.528731198194765,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
      },
      {
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        2.528731198194765,
        0.1470905524127557,
        -0.1470905524127557,
      },
      {
        -0.008416385185447427,
        0.07051733851899882,
        0.1149203802569082,
        0.0000000000000000,
        0.1470905524127557,
        0.760031883566609,
        -0.003979736727037247,
      },
      {
        -0.008416385185447427,
        0.07051733851899882,
        0.1149203802569082,
        0.0000000000000000,
        -0.1470905524127557,
        -0.003979736727037247,
        0.760031883566609,
      },
    };
}
} // namespace test

TEST_CASE("Kinetic") {
    using test_pt = simde::aos_t_e_aos;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("Kinetic"));

    // Get basis set
    auto mol  = test::water_molecule();
    auto aobs = test::water_sto3g_basis_set();

    // Make AOS object
    simde::type::aos aos(aobs);

    // Make Operator
    simde::type::t_e_type op{};

    // Make BraKet Input
    chemist::braket::BraKet braket(aos, op, aos);

    // Call module
    auto T = mm.at("Kinetic").run_as<test_pt>(braket);
    auto corr = test::correct_value();

#define DOWNCAST static_cast<const tensorwrapper::buffer::Eigen<double, 2>&>

    auto t = DOWNCAST(T.buffer());
    std::cout << t.value() << std::endl;

    auto t_corr = DOWNCAST(corr.buffer());
    std::cout << t_corr.value() << std::endl;

#undef DOWNCAST

    // Check output
    REQUIRE(T == corr);
}
