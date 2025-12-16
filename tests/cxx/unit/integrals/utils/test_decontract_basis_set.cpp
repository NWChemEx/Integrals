/*
 * Copyright 2025 NWChemEx-Project
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
#include <integrals/property_types.hpp>

TEST_CASE("DecontractBasisSet") {
    using test_pt = integrals::property_types::decontract_basis_set;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    REQUIRE(mm.count("Decontract Basis Set"));

    // Get basis set
    auto aobs = test::water_sto3g_basis_set();

    // Get module
    auto& mod = mm.at("Decontract Basis Set");

    // Call module
    auto decontracted_aobs = mod.run_as<test_pt>(aobs);

    // Check output
    auto corr = test::water_decontracted_sto3g_basis_set();
    REQUIRE(decontracted_aobs == corr);
}
