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

#include "../utils/rank2_shell_norm.hpp"
#include "libint.hpp"
#include <cmath>
#include <integrals/property_types.hpp>

namespace integrals::libint {
namespace {

const auto desc = "";

} // namespace

using decontract_pt = integrals::property_types::DecontractBasisSet;
using eri4_pt       = simde::ERI4;
using pt            = integrals::property_types::PrimitivePairEstimator;

MODULE_CTOR(CauchySchwarzPrimitiveEstimator) {
    satisfies_property_type<pt>();
    description(desc);
    // TODO: Add citation for Chemist paper
    add_submodule<decontract_pt>("Decontract Basis Set");
    add_submodule<eri4_pt>("ERI4");
}

MODULE_RUN(CauchySchwarzPrimitiveEstimator) {
    const auto&& [bra_basis, ket_basis] = pt::unwrap_inputs(inputs);

    auto& to_prims_mod    = submods.at("Decontract Basis Set");
    const auto& bra_prims = to_prims_mod.run_as<decontract_pt>(bra_basis);
    const auto& ket_prims = to_prims_mod.run_as<decontract_pt>(ket_basis);
    const auto n_bra      = bra_prims.n_primitives();
    const auto n_ket      = ket_prims.n_primitives();

    // Should always be true, but we check for sanity
    assert(n_bra == bra_basis.n_primitives());
    assert(n_ket == ket_basis.n_primitives());

    // TODO: Only get the hyper diagonal we actually need
    simde::type::aos_squared bra(bra_prims, ket_prims);
    simde::type::v_ee_type v_ee{};
    chemist::braket::BraKet mnls(bra, v_ee, bra);
    const auto& prim4 = submods.at("ERI4").run_as<eri4_pt>(mnls);

    using tensorwrapper::buffer::make_contiguous;
    const auto& eris = make_contiguous(prim4.buffer());

    // TODO: Use floating point type of the basis sets
    using float_type = double;
    std::vector<float_type> buffer(n_bra * n_ket, 0.0);

    using iter_type = std::decay_t<decltype(n_bra)>; // Type of the loop index

    // We loop over the decontracted basis, but make sure to grab the
    // contraction coefficients from the real basis set. By looping over the
    // decontracted basis we can make this look like the "usual" shell pair
    // loop
    // N.b. there's one shell per primitive in the real basis
    // N.b. since the basis is decontracted the number of Cartesian/spherical
    // AOs in a shell is the same as the number of Cartesian/spherical
    // primitives in that shell

    std::array<iter_type, 4> offsets{0, 0, 0, 0};
    std::array<iter_type, 4> naos{0, 0, 0, 0};
    for(iter_type shell_i = 0; shell_i < bra_prims.n_shells(); ++shell_i) {
        const auto bra_prim_i = bra_basis.primitive(shell_i);
        const auto bra_coeff  = std::fabs(bra_prim_i.coefficient());
        naos[0]               = bra_prims.shell(shell_i).size();
        naos[2]               = naos[0];

        const auto bra_shell_offset = shell_i * n_ket;
        offsets[1]                  = 0;
        offsets[3]                  = 0;
        for(iter_type shell_j = 0; shell_j < ket_prims.n_shells(); ++shell_j) {
            const auto ket_prim_i = ket_basis.primitive(shell_j);
            const auto ket_coeff  = std::fabs(ket_prim_i.coefficient());
            naos[1]               = ket_prims.shell(shell_j).size();
            naos[3]               = naos[1];

            auto C_ab       = bra_coeff * ket_coeff;
            auto shell_norm = utils::rank2_shell_norm(eris, offsets, naos);
            buffer[bra_shell_offset + shell_j] = C_ab * std::sqrt(shell_norm);

            offsets[1] += naos[1];
            offsets[3] += naos[3];
        }
        offsets[0] += naos[0];
        offsets[2] += naos[2];
    }

    tensorwrapper::shape::Smooth shape({n_bra, n_ket});
    tensorwrapper::buffer::Contiguous tw_buffer(std::move(buffer), shape);
    simde::type::tensor rv(shape, std::move(tw_buffer));

    auto result = results();
    return pt::wrap_results(result, rv);
}

} // namespace integrals::libint