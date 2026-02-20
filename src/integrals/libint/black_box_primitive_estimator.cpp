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

#include "libint.hpp"
#include <cmath>
#include <integrals/property_types.hpp>

namespace integrals::libint {
namespace {

const auto desc = "";

}

using pt = integrals::property_types::PrimitivePairEstimator;

MODULE_CTOR(BlackBoxPrimitiveEstimator) {
    satisfies_property_type<pt>();
    description(desc);
    // TODO: Add citation for Chemist paper
}

MODULE_RUN(BlackBoxPrimitiveEstimator) {
    const auto&& [bra_basis, ket_basis] = pt::unwrap_inputs(inputs);

    const auto n_shells_bra = bra_basis.n_shells();
    const auto n_prims_bra  = bra_basis.n_primitives();
    const auto n_shells_ket = ket_basis.n_shells();
    const auto n_prims_ket  = ket_basis.n_primitives();

    // TODO: Our basis really needs to handle normalization better...
    auto bra = make_libint_basis_set(bra_basis);
    auto ket = make_libint_basis_set(ket_basis);

    // TODO: Use floating point type of the basis sets
    using float_type = double;
    std::vector<float_type> buffer(n_prims_bra * n_prims_ket, 0.0);

    using iter_type = std::decay_t<decltype(n_bra)>; // Type of the loop index
    iter_type bra_counter = 0;

    for(iter_type bra_shell_i = 0; bra_shell_i < n_shells_bra; ++bra_shell_i) {
        const auto& bra_shell = bra.at(bra_shell_i);
        assert(bra_shell.contr.size() == 1);
        const auto& bra_coeff        = bra_shell.contr[0].coeff;
        const auto& bra_alpha        = bra_shell.contr[0].alpha;
        const auto n_prims_bra_shell = bra_coeff.size();

        for(iter_type bra_prim_i = 0; bra_prim_i < n_prims_bra_shell;
            ++bra_prim_i) {
            const auto bra_zeta   = bra_alpha[bra_prim_i];
            const auto bra_coeff  = std::fabs(bra_coeff[bra_prim_i]);
            const auto bra_center = bra_shell.O;
            const auto bra_offset = (bra_counter + bra_prim_i) * n_prims_ket;

            iter_type ket_counter = 0;
            for(iter_type ket_shell_i = 0; ket_shell_i < n_shells_ket;
                ++ket_shell_i) {
                const auto& ket_shell = ket.at(ket_shell_i);
                assert(ket_shell.contr.size() == 1);
                const auto& ket_coeff        = ket_shell.contr[0].coeff;
                const auto& ket_alpha        = ket_shell.contr[0].alpha;
                const auto n_prims_ket_shell = ket_coeff.size();

                for(iter_type ket_prim_i = 0; ket_prim_i < n_prims_ket_shell;
                    ++ket_prim_i) {
                    const auto ket_zeta   = ket_alpha[ket_prim_i];
                    const auto ket_coeff  = std::fabs(ket_coeff[ket_prim_i]);
                    const auto ket_center = ket_shell.O;
                    const auto ket_offset = ket_counter + ket_prim_i;

                    // This is "K bar" in Eq. 11 in the SI of the Chemist paper
                    const auto dx  = bra_center.x - ket_center.x;
                    const auto dy  = bra_center.y - ket_center.y;
                    const auto dz  = bra_center.z - ket_center.z;
                    const auto dr2 = dx * dx + dy * dy + dz * dz;

                    const auto ratio =
                      bra_zeta * ket_zeta / (bra_zeta + ket_zeta);
                    const auto coeff = bra_coeff * ket_coeff;
                    buffer[bra_offset + ket_offset] =
                      coeff * std::exp(-1.0 * ratio * dist2);
                }
                ket_counter += n_prims_ket_shell;
            }
            bra_counter += n_prims_bra_shell;
        }
        tensorwrapper::shape::Smooth shape({n_bra, n_ket});
        tensorwrapper::buffer::Contiguous tw_buffer(std::move(buffer), shape);
        simde::type::tensor rv(shape, std::move(tw_buffer));

        auto result = results();
        return pt::wrap_results(result, rv);
    }

} // namespace integrals::libint