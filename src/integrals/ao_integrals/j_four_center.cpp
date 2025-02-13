/*
 * Copyright 2024 NWChemEx-Project
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

using pt    = simde::aos_j_e_aos;
using pt_4c = simde::ERI4;

namespace {

auto desc = R"(
Four-Center J Builder
---------------------
)";

}
MODULE_CTOR(JFourCenter) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<pt_4c>("Four-center ERI");
}

MODULE_RUN(JFourCenter) {
    const auto&& [braket] = pt::unwrap_inputs(inputs);
    // TODO: avoid copying AOs
    simde::type::aos bra_e0 = braket.bra();
    const auto& j_hat       = braket.op();
    simde::type::aos ket_e0 = braket.ket();
    const auto& rho         = j_hat.rhs_particle();
    simde::type::aos aos_e1 = rho.basis_set();
    const auto& p           = rho.value();
    auto& eri_mod           = submods.at("Four-center ERI");

    // auto aos2_v_aos2 = (bra_e0 * ket_e0 | v_ee | aos_e1 * aos_e1);
    simde::type::v_ee_type v_ee;
    simde::type::aos_squared e0_pair(bra_e0, ket_e0);
    simde::type::aos_squared e1_pair(aos_e1, aos_e1);
    chemist::braket::BraKet aos2_v_aos2(e0_pair, v_ee, e1_pair);
    const auto& I = eri_mod.run_as<pt_4c>(std::move(aos2_v_aos2));

    simde::type::tensor j;
    j.multiplication_assignment("i,j", p("k,l"), I("i,j,k,l"));

    auto rv = results();
    return pt::wrap_results(rv, std::move(j));
}

} // namespace integrals::ao_integrals