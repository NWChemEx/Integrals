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

#include "detail_/make_libint_basis_set.hpp"
#include "libint.hpp"
#include <cmath>
#include <integrals/property_types.hpp>

namespace integrals::libint {
namespace {

const auto desc = R"(
Libint Black Box Primitive Pair Estimator 
=========================================

This module computes the matrix : math :`K_{ij}` where:

.. math:: 

   K_{ij} = c_i c_j \exp\left(-\frac{\zeta_i \zeta_j}{\zeta_i + \zeta_j} 
   |\mathbf{R}_i - \mathbf{R}_j|^2\right)

This is how Libint2 estimates the contribution of a pair of primitives to
an integral.

N.B. The algorithm assumes that the bra and ket are different. If they are the
same, we can save time by using the fact that the matrix is symmetric.

)";

// Computes square of the distance between points a and b
// (T should be libint::Atom like)
template<typename T>
auto distance_squared(T&& a, T&& b) {
    auto dx = a[0] - b[0];
    auto dy = a[1] - b[1];
    auto dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
}

template<typename T>
auto compute_k(T zeta_i, T zeta_j, T coeff_i, T coeff_j, T dr2) {
    const auto num   = -zeta_i * zeta_j;
    const auto denom = zeta_i + zeta_j;
    const auto ratio = num / denom;
    return coeff_i * coeff_j * std::exp(ratio * dr2);
}

} // namespace
using pt = integrals::property_types::PrimitivePairEstimator;

MODULE_CTOR(BlackBoxPrimitiveEstimator) {
    satisfies_property_type<pt>();
    description(desc);
    // TODO: Add citation for Chemist paper
}

MODULE_RUN(BlackBoxPrimitiveEstimator) {
    const auto&& [bra, ket] = pt::unwrap_inputs(inputs);

    using iter_type    = std::size_t;              // Type of each loop index
    using index_array  = std::array<iter_type, 2>; // Type of a pair of indices
    using index_vector = std::vector<iter_type>; // Type of a vector of indices
    using float_type   = double;                 // TODO: Get from basis sets

    index_array shells{0, 0}; // shells[0]/shells[1] indexes bra/ket shell
    index_array n_shells{bra.n_shells(), ket.n_shells()}; // Number of shells

    index_array prims{0, 0}; // prims[0]/prims[1] indexes bra/ket primitive
    index_array n_prims{bra.n_primitives(), ket.n_primitives()};
    index_array offsets{0, 0};   // Offset to the first primitive of the shell
    index_vector abs_prim{0, 0}; // Absolute indices of the primitives

    // Will be the result
    tensorwrapper::shape::Smooth shape({n_prims[0], n_prims[1]});
    std::vector<float_type> data(shape.size(), 0.0);
    tensorwrapper::buffer::Contiguous buffer(std::move(data), shape);

    // For now use the libint basis sets because they're properly normalized
    // TODO: Our basis really needs to handle normalization better...
    auto bra_libint = detail_::make_libint_basis_set(bra);
    auto ket_libint = detail_::make_libint_basis_set(ket);

    for(shells[0] = 0; shells[0] < n_shells[0]; ++shells[0]) {
        const auto& bra_shell = bra_libint.at(shells[0]);
        assert(bra_shell.contr.size() == 1); // No general contraction support
        const auto& bra_coeff        = bra_shell.contr[0].coeff;
        const auto& bra_alpha        = bra_shell.alpha;
        const auto n_prims_bra_shell = bra_coeff.size();

        for(prims[0] = 0; prims[0] < n_prims_bra_shell; ++prims[0]) {
            const auto zeta0      = bra_alpha[prims[0]];
            const auto coeff0     = std::fabs(bra_coeff[prims[0]]);
            const auto bra_center = bra_shell.O;
            abs_prim[0]           = offsets[0] + prims[0];

            offsets[1] = 0;
            for(shells[1] = 0; shells[1] < n_shells[1]; ++shells[1]) {
                const auto& ket_shell = ket_libint.at(shells[1]);
                assert(ket_shell.contr.size() == 1); // No general contractions
                const auto& ket_coeff        = ket_shell.contr[0].coeff;
                const auto& ket_alpha        = ket_shell.alpha;
                const auto n_prims_ket_shell = ket_coeff.size();

                for(prims[1] = 0; prims[1] < n_prims_ket_shell; ++prims[1]) {
                    const auto zeta1      = ket_alpha[prims[1]];
                    const auto coeff1     = std::fabs(ket_coeff[prims[1]]);
                    const auto ket_center = ket_shell.O;
                    abs_prim[1]           = offsets[1] + prims[1];

                    // This is "K bar" in Eq. 11 in the SI of the Chemist paper
                    auto dr2 = distance_squared(bra_center, ket_center);
                    auto k01 = compute_k(zeta0, zeta1, coeff0, coeff1, dr2);
                    buffer.set_elem(abs_prim, k01);
                } // loop over ket primitives

                offsets[1] += n_prims_ket_shell;
            } // loop over ket shells

            // We ran over all ket primitives, so counter should be done too
            assert(offsets[1] == n_prims[1]);
        } // loop over bra primitives

        offsets[0] += n_prims_bra_shell;
    } // loop over bra shells

    // We ran over all bra primitives, so counter should be done too
    assert(offsets[0] == n_prims[0]);

    simde::type::tensor rv(shape, std::move(buffer));

    auto result = results();
    return pt::wrap_results(result, rv);
}

} // namespace integrals::libint