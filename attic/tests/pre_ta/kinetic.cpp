/*
 * Copyright 2023 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
#include <integrals/libint_integral.hpp>
#include <property_types/aointegral.hpp>

using namespace integrals::libint;

// Computes the kinetic energy integrals for water in STO-3G
static BlockTensor corr{
  {{
     0,
     0,
   },
   {
     29.0031999455395848, -0.1680109393164922, 0.0000000000000000,
     0.0000000000000000,  0.0000000000000000,  -0.0084163851854474,
     -0.0084163851854474, -0.1680109393164923, 0.8081279549303477,
     -0.0000000000000000, 0.0000000000000000,  0.0000000000000000,
     0.0705173385189988,  0.0705173385189988,  0.0000000000000000,
     -0.0000000000000000, 2.5287311981947651,  0.0000000000000000,
     0.0000000000000000,  0.1149203802569082,  0.1149203802569082,
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
     2.5287311981947642,  0.0000000000000000,  0.0000000000000000,
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
     0.0000000000000000,  0.0000000000000000,  2.5287311981947651,
     0.1470905524127557,  -0.1470905524127557, -0.0084163851854474,
     0.0705173385189988,  0.1149203802569082,  0.0000000000000000,
     0.1470905524127557,  0.7600318835666090,  -0.0039797367270372,
     -0.0084163851854474, 0.0705173385189988,  0.1149203802569082,
     0.0000000000000000,  -0.1470905524127557, -0.0039797367270372,
     0.7600318835666090,
   }}};

TEST_CASE("Testing Libint's Kinetic Energy Integrals class") {
    auto[molecule, bs]                          = make_molecule();
    std::array<chemist::AOBasisSet, 2> bases = {bs, bs};
    using integral_type = property_types::AOIntegral<2, double>;
    sde::ModuleManager mm;
    SECTION("Direct implementation") {
        load_modules(mm);
        auto[Ints] = mm.at("Kinetic").run_as<integral_type>(molecule, bases,
                                                            std::size_t{0});
        compare_integrals(Ints, corr);
    }
    SECTION("Core implementation") {
        mm.add_module("KineticCore", std::make_shared<Kinetic>(
                                       detail_::implementation_type::core));
        auto[Ints] = mm.at("KineticCore")
                       .run_as<integral_type>(molecule, bases, std::size_t{0});
        compare_integrals(Ints, corr);
    }
}
