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
#include "test_error.hpp"
#include <integrals/property_types.hpp>

using pt = integrals::property_types::PrimitivePairEstimator;
using namespace integrals::libint::test;
TEST_CASE("CauchySchwarzPrimitiveEstimator") {
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    integrals::set_defaults(mm);
    auto& mod                     = mm.at("CauchySchwarz Estimator");
    auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

    auto Q_ab = mod.run_as<pt>(bra0, bra1);
    auto Q_cd = mod.run_as<pt>(ket0, ket1);
    // std::cout << Q_ab << std::endl;
}