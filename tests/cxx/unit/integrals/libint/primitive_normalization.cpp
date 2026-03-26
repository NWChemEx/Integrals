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
#include <integrals/property_types.hpp>

using namespace integrals::testing;

using pt = integrals::property_types::Normalize<simde::type::ao_basis_set>;

namespace {

/* Renormalized contraction coefficients for water STO-3G, matching libint's
 * renorm() convention.
 *
 * For each contracted shell, libint applies:
 *   1. Per-primitive factor N_p = sqrt(2^l (2*zeta_p)^(l+3/2) /
 *      (sqrt(pi^3) * (2l-1)!!))
 *   2. Contracted-shell unit-norm scaling: divide all d_p*N_p by
 *      sqrt(<phi|phi>) where <phi|phi> is the contracted shell self-overlap.
 *
 * The returned values are |d_p * N_p / sqrt(<phi|phi>)| for each primitive,
 * replicated once per AO component within that primitive.
 *
 * Shell order matches the contracted basis (not decontracted):
 *   O s1 : 3 primitives (zeta = 130.7093200, 23.8088610, 6.4436083)
 *   O s2 : 3 primitives (zeta =   5.0331513,  1.1695961, 0.3803890)
 *   O p  : 3 primitives × 3 AOs (zeta = 5.0331513, 1.1695961, 0.3803890)
 *   H1 s : 3 primitives (zeta = 3.425250914, 0.6239137298, 0.1688554040)
 *   H2 s : 3 primitives (same as H1)
 */
std::vector<double> corr_water_norms() {
    return {
      // O s1
      4.251943282943720,
      4.112293718431184,
      1.281622532581341,
      // O s2
      0.239413002994477,
      0.320234229133891,
      0.241685570753216,
      // O p (3 primitives × 3 AOs = 9 entries, same value per AO)
      1.675450118114190,
      1.675450118114190,
      1.675450118114190,
      1.053568007994846,
      1.053568007994846,
      1.053568007994846,
      0.166902898075749,
      0.166902898075749,
      0.166902898075749,
      // H1 s
      0.276934355079052,
      0.267838851609479,
      0.083473671129841,
      // H2 s
      0.276934355079052,
      0.267838851609479,
      0.083473671129841,
    };
}

} // namespace

TEST_CASE("PrimitiveNormalization") {
    auto mm   = initialize_integrals();
    auto& mod = mm.at("Primitive Normalization");

    SECTION("H2O/STO-3G") {
        auto aobs  = water_sto3g_basis_set();
        auto norms = mod.run_as<pt>(aobs);
        auto corr  = corr_water_norms();

        REQUIRE(norms.size() == corr.size());
        for(std::size_t i = 0; i < norms.size(); ++i) {
            REQUIRE(norms[i] == Catch::Approx(corr[i]).epsilon(1.0e-10));
        }
    }
}
