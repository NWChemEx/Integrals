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

/* Reference normalization factors 1/sqrt((p|p)) for water STO-3G, computed
 * analytically from the Gaussian overlap formula.
 *
 * The decontracted basis order follows DecontractBasisSet:
 *   O s1 : 3 primitives (zeta = 130.7093200, 23.8088610, 6.4436083)
 *   O s2 : 3 primitives (zeta =   5.0331513,  1.1695961, 0.3803890)
 *   O p  : 3 primitives × 3 AOs (zeta = 5.0331513, 1.1695961, 0.3803890)
 *   H1 s : 3 primitives (zeta = 3.425250914, 0.6239137298, 0.1688554040)
 *   H2 s : 3 primitives (same as H1)
 *
 * For s-type: (p|p) = (pi/(2*zeta))^(3/2)
 * For p-type (e.g. x*exp(-zeta*r^2)):
 *   (p|p) = integral x^2 exp(-2*zeta*r^2) d^3r = (1/(4*zeta)) *
 * (pi/(2*zeta))^(3/2) Normalization factor = 1/sqrt((p|p))
 */
std::vector<double> corr_water_norms() {
    return {
      // O s1
      27.551167600757328,
      7.681818767185153,
      2.882417868807197,
      // O s2
      2.394914875721071,
      0.801561825779394,
      0.345208161163625,
      // O p (3 primitives × 3 AOs = 9 entries, same value per AO)
      10.745832583524919,
      10.745832583524919,
      10.745832583524919,
      1.733744024405248,
      1.733744024405248,
      1.733744024405248,
      0.425818989418287,
      0.425818989418287,
      0.425818989418287,
      // H1 s
      1.794441833790094,
      0.500326492211116,
      0.187735461846361,
      // H2 s
      1.794441833790094,
      0.500326492211116,
      0.187735461846361,
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
