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

#include "test_common.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/aointegral.hpp>

using namespace integrals::libint;

// Computes the nuclear-electron energy integrals for water in STO-3G
static BlockTensor corr{
  {{
     0,
     0,
   },
   {
     -61.5805952694322016, -7.4108218563311627, -0.0144738837457361,
     0.0000000000000000,   0.0000000000000000,  -1.2316855721424882,
     -1.2316855721424882,  -7.4108218563311627, -10.0090711420700309,
     -0.1768908347336432,  0.0000000000000000,  0.0000000000000000,
     -2.9772268535781348,  -2.9772268535781343, -0.0144738837457361,
     -0.1768908347336432,  -9.9440433416987624, 0.0000000000000000,
     0.0000000000000000,   -1.4717933387129607, -1.4717933387129607,
     0.0000000000000000,   0.0000000000000000,  0.0000000000000000,
     -9.8758759950909436,  0.0000000000000000,  0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -9.9875499350885555,
     -1.8222369134761305,  1.8222369134761305,  -1.2316855721424882,
     -2.9772268535781339,  -1.4717933387129609, 0.0000000000000000,
     -1.8222369134761305,  -5.3002032522950211, -1.0671710804724368,
     -1.2316855721424882,  -2.9772268535781343, -1.4717933387129609,
     0.0000000000000000,   1.8222369134761305,  -1.0671710804724368,
     -5.3002032522950193,
   }}};

TEST_CASE("Testing Libint Nuclear-Electron Attraction Integrals") {
    using integral_type = property_types::AOIntegral<2, double>;
    sde::ModuleManager mm;
    load_modules(mm);
    auto[molecule, bs]                          = make_molecule();
    std::array<chemist::AOBasisSet, 2> bases = {bs, bs};
    auto[Ints] =
      mm.at("Nuclear").run_as<integral_type>(molecule, bases, std::size_t{0});
    compare_integrals(Ints, corr);
}
