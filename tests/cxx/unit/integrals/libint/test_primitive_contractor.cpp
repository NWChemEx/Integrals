/*
 * Copyright 2026 NWChemEx-Project
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
using tensorwrapper::operations::approximately_equal;

TEST_CASE("Primitive Contractor ERI4") {
    using test_pt = simde::ERI4;

    auto mm = initialize_integrals();
    REQUIRE(mm.count("Primitive Contractor ERI4"));
    REQUIRE(mm.count("ERI4"));

    auto aobs = water_sto3g_basis_set();
    for(const auto& type :
        {chemist::ShellType::cartesian, chemist::ShellType::pure}) {
        for(const auto& l : {1, 2, 3}) {
            aobs.shell(2).l()    = l;
            aobs.shell(2).pure() = type;
            simde::type::aos aos(aobs);
            simde::type::aos_squared aos_squared(aos, aos);

            simde::type::v_ee_type op{};

            chemist::braket::BraKet braket(aos_squared, op, aos_squared);
            auto corr_mod = mm.at("ERI4");
            auto test_mod = mm.at("Primitive Contractor ERI4");

            auto T_contracted = test_mod.run_as<test_pt>(braket);
            auto T_eri4       = corr_mod.run_as<test_pt>(braket);

            REQUIRE(approximately_equal(T_contracted, T_eri4, 1.0e-14));
        }
    }
}
