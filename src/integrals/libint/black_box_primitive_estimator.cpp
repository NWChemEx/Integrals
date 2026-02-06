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

    const auto n_bra = bra_basis.n_primitives();
    const auto n_ket = ket_basis.n_primitives();

    // TODO: Use floating point type of the basis sets
    using float_type = double;
    std::vector<float_type> buffer(n_bra * n_ket, 0.0);

    using iter_type = std::decay_t<decltype(n_bra)>; // Type of the loop index

    for(iter_type bra_i = 0; bra_i < n_bra; ++bra_i) {
        const auto bra_prim_i = bra_basis.primitive(bra_i);
        const auto bra_zeta   = bra_prim_i.exponent();
        const auto bra_coeff  = std::fabs(bra_prim_i.coefficient());
        const auto bra_center = bra_prim_i.center();
        const auto bra_offset = bra_i * n_ket;

        for(iter_type ket_i = 0; ket_i < n_ket; ++ket_i) {
            const auto ket_prim_i = ket_basis.primitive(ket_i);
            const auto ket_zeta   = ket_prim_i.exponent();
            const auto ket_coeff  = std::fabs(ket_prim_i.coefficient());
            const auto ket_center = ket_prim_i.center();

            // This is "K bar" in Eq. 11 in the SI of the Chemist paper
            const auto dist  = (bra_center - ket_center).magnitude();
            const auto dist2 = dist * dist;
            const auto ratio = bra_zeta * ket_zeta / (bra_zeta + ket_zeta);
            const auto coeff = bra_coeff * ket_coeff;
            buffer[bra_offset + ket_i] = coeff * std::exp(-1.0 * ratio * dist2);
        }
    }

    tensorwrapper::shape::Smooth shape({n_bra, n_ket});
    tensorwrapper::buffer::Contiguous tw_buffer(std::move(buffer), shape);
    simde::type::tensor rv(shape, std::move(tw_buffer));

    auto result = results();
    return pt::wrap_results(result, rv);
}

} // namespace integrals::libint