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

    auto corr_mod = mm.at("ERI4");
    auto test_mod = mm.at("Primitive Contractor ERI4");

    auto aobs               = water_sto3g_basis_set();
    const auto is_cartesian = chemist::ShellType::cartesian;
    const auto is_pure      = chemist::ShellType::pure;
    for(const auto& type : {is_cartesian, is_pure}) {
        const auto type_str =
          (type == is_cartesian) ? "Cartesian" : "Spherical";
        SECTION(type_str) {
            for(const auto& l : {1, 2, 3}) {
                SECTION("Angular Momentum: " + std::to_string(l)) {
                    aobs.shell(2).l()    = l;
                    aobs.shell(2).pure() = type;
                    simde::type::aos aos(aobs);
                    simde::type::aos_squared aos_squared(aos, aos);
                    simde::type::v_ee_type op{};
                    chemist::braket::BraKet braket(aos_squared, op,
                                                   aos_squared);

                    for(const auto& thresh : {1e-16, 1e-8, 1e-6}) {
                        SECTION("Screening Threshold: " +
                                std::to_string(thresh)) {
                            auto test_mod_copy = test_mod.unlocked_copy();
                            test_mod_copy.change_input("Screening Threshold",
                                                       thresh);

                            auto corr_mod_copy = corr_mod.unlocked_copy();
                            corr_mod_copy.change_input("Threshold", thresh);
                            auto T_contracted =
                              test_mod_copy.run_as<test_pt>(braket);
                            auto T_eri4 = corr_mod_copy.run_as<test_pt>(braket);

                            REQUIRE(approximately_equal(T_contracted, T_eri4,
                                                        1.0e-14));
                        }
                    }
                }
            }
        }
    }
}
