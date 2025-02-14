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

#include "ao_integrals.hpp"

namespace integrals::ao_integrals {

using pt    = simde::ERI3;
using pt_2c = simde::ERI2;

namespace {

auto desc = R"(
Three-index ERI with Coulomb metric transformation
---------------------
)";

}

MODULE_CTOR(DFIntegral) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<pt>("Three-center ERI");
    add_submodule<pt_2c>("Coulomb Metric");
}

MODULE_RUN(DFIntegral) {
    const auto& [braket] = pt::unwrap_inputs(inputs);
    auto bra             = braket.bra();
    auto ket             = braket.ket();
    auto& op             = braket.op();
    auto& eri2_mod       = submods.at("Coulomb Metric");
    auto& eri3_mod       = submods.at("Three-center ERI");

    chemist::braket::BraKet aux_v_aux(bra, op, bra);
    const auto& M = eri2_mod.run_as<pt_2c>(aux_v_aux);
    const auto& I = eri3_mod.run_as<pt>(braket);

    // Failing at the moment
    simde::type::tensor L;
    // L.multiplication_assignment("i,k,l", M("i,j"), I("j,k,l"));

    auto rv = results();
    return pt::wrap_results(rv, std::move(L));
}

} // namespace integrals::ao_integrals