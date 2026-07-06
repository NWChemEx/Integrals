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
using pt_3c = simde::ERI3;
using aos_t = simde::type::aos;

namespace {

auto desc = R"(
Density Fitted J Builder
---------------------
)";

}
MODULE_CTOR(JDensityFitted) {
    description(desc);
    satisfies_property_type<pt>();
    add_submodule<pt_3c>("DF ERI");
    add_input<aos_t>("Auxiliary Basis Set");
}

MODULE_RUN(JDensityFitted) {
    const auto&& [braket] = pt::unwrap_inputs(inputs);

    auto bra          = braket.bra();
    auto ket          = braket.ket();
    const auto& j_hat = braket.op();
    const auto& rho   = j_hat.get_rhs_particle();
    const auto& p     = rho.value();
    auto rho_aos      = rho.basis_set();
    auto aux          = inputs.at("Auxiliary Basis Set").value<aos_t>();
    auto& eri_mod     = submods.at("DF ERI");

    simde::type::v_ee_type v_ee;
    simde::type::aos_squared ij_pair(bra, ket);
    simde::type::aos_squared kl_pair(rho_aos, rho_aos);
    chemist::braket::BraKet aux_v_ij(aux, v_ee, ij_pair);
    chemist::braket::BraKet aux_v_kl(aux, v_ee, kl_pair);
    const auto& I_akl = eri_mod.run_as<pt_3c>(std::move(aux_v_kl));
    const auto& I_aij = eri_mod.run_as<pt_3c>(std::move(aux_v_ij));

    simde::type::tensor j;
    j("a")   = p("k,l") * I_akl("a,k,l");
    j("i,j") = j("a") * I_aij("a,i,j");

    auto rv = results();
    return pt::wrap_results(rv, std::move(j));
}

} // namespace integrals::ao_integrals
