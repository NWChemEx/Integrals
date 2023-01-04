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

#include "H2O_STO3G_STG[1].hpp"
#include "H2O_STO3G_Yukawa[1].hpp"
#include "test_common.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/aointegral.hpp>

using namespace integrals::libint;

// Computes the ERI integrals for water in STO-3G
TEST_CASE("Testing Libint's integrals over Slater-type geminals") {
    using integral_type = property_types::AOIntegral<4, double>;
    sde::ModuleManager mm;
    load_modules(mm);
    auto[molecule, bs]                          = make_molecule();
    std::array<chemist::AOBasisSet, 4> bases = {bs, bs, bs, bs};
    auto[STG4_Ints] = mm.at("STG4").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(STG4_Ints, stg1ref, 0.0, 1e-12);
    auto[Yukawa4_Ints] = mm.at("Yukawa4").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Yukawa4_Ints, yukawa1ref, 0.0, 1e-12);
}
