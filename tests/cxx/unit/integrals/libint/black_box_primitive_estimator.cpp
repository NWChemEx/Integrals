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
#include "testing.hpp"
#include <integrals/property_types.hpp>

using pt = integrals::property_types::PrimitivePairEstimator;

namespace {

std::vector<double> h2_h2_corr{
  0.023817430147884473,   0.08261663936778646, 0.06861998972363427,
  1.4555047643750448e-05, 0.00844595180076824, 0.03423471319097009,
  0.08261663936778646,    0.28657621993836907, 0.23802538347833696,
  0.00844595180076824,    0.07444365105885908, 0.13404298325096972,
  0.06861998972363427,    0.23802538347833696, 0.19769987611740356,
  0.03423471319097009,    0.13404298325096972, 0.1372685140907819,
  1.4555047643750448e-05, 0.00844595180076824, 0.03423471319097009,
  0.023817430147884473,   0.08261663936778646, 0.06861998972363427,
  0.00844595180076824,    0.07444365105885908, 0.13404298325096972,
  0.08261663936778646,    0.28657621993836907, 0.23802538347833696,
  0.03423471319097009,    0.13404298325096972, 0.1372685140907819,
  0.06861998972363427,    0.23802538347833696, 0.19769987611740356};

} // namespace

TEST_CASE("BlackBoxPrimitiveEstimator") {
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);
    auto& mod     = mm.at("Black Box Primitive Pair Estimator");
    auto h2_basis = test::h2_sto3g_basis_set();

    SECTION("H2 STO-3G") {
        const auto& Kij = mod.run_as<pt>(h2_basis, h2_basis);
        tensorwrapper::shape::Smooth shape{6, 6};
        tensorwrapper::buffer::Contiguous corr(h2_h2_corr, shape);
        REQUIRE(Kij.buffer().approximately_equal(corr, 1e-8));
    }

    // Tests bra != ket
    SECTION("Block of H2 STO-3G") {
        simde::type::ao_basis_set h0_basis;
        h0_basis.add_center(h2_basis.at(0));
        const auto& Kij = mod.run_as<pt>(h0_basis, h2_basis);
        tensorwrapper::shape::Smooth shape{3, 6};
        std::vector<double> buffer(h2_h2_corr.begin(), h2_h2_corr.begin() + 18);
        tensorwrapper::buffer::Contiguous corr(buffer, shape);
        REQUIRE(Kij.buffer().approximately_equal(corr, 1e-8));
    }
}