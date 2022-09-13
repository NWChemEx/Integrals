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

#include "H2O_STO3G_DF.hpp"
#include "test_common.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/aointegral.hpp>

using namespace integrals::libint;

// Computes the three-center, two-electron integrals for water in sto-3g
// Note: as with the metric tensor using sto-3g instead of a proper fitting
// basis  should be fine for testing purposes.

TEST_CASE("Testing DF3C2E") {
    using integral_type = property_types::AOIntegral<3, double>;
    sde::ModuleManager mm;
    load_modules(mm);
    auto[molecule, bs]                          = make_molecule();
    std::array<chemist::AOBasisSet, 3> bases = {bs, bs, bs};
    auto[Ints] =
      mm.at("ERI3").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}
