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
    auto bra              = braket.bra();
    auto ket              = braket.ket();
    const auto& j_hat     = braket.op();
    const auto& rho       = j_hat.rhs_particle();
    const auto& p         = rho.value();
    auto rho_aos          = rho.basis_set();
    auto& eri_mod         = submods.at("Four-center ERI");

    simde::type::v_ee_type v_ee;
    simde::type::aos_squared ij_pair(bra, ket);
    simde::type::aos_squared kl_pair(rho_aos, rho_aos);
    chemist::braket::BraKet ij_v_kl(ij_pair, v_ee, kl_pair);
    const auto& I = eri_mod.run_as<pt_4c>(std::move(ij_v_kl));

    simde::type::tensor j;
    j("i,j") = p("k,l") * I("i,j,k,l");

    auto rv = results();
    return pt::wrap_results(rv, std::move(j));
}

} // namespace integrals::ao_integrals
