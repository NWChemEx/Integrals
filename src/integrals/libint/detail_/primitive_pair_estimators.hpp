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

#pragma once
#include "make_libint_basis_set.hpp"
#include <simde/simde.hpp>
#include <vector>
namespace integrals::libint::detail_ {

inline auto k_ij(const simde::type::ao_basis_set& basis0,
                 const simde::type::ao_basis_set& basis1) {
    auto distance_squared = [](auto&& a, auto&& b) {
        auto dx = a.x() - b.x();
        auto dy = a.y() - b.y();
        auto dz = a.z() - b.z();
        return dx * dx + dy * dy + dz * dz;
    };

    auto nprims0   = basis0.n_primitives();
    auto nprims1   = basis1.n_primitives();
    using vector_t = std::vector<double>;
    using matrix_t = std::vector<vector_t>;
    matrix_t K(nprims0, vector_t(nprims1, 0.0));
    for(std::size_t i = 0; i < nprims0; ++i) {
        auto alpha0 = basis0.primitive(i).exponent();
        auto r0     = basis0.primitive(i).center();
        for(std::size_t j = 0; j < nprims1; ++j) {
            auto alpha1 = basis1.primitive(j).exponent();
            auto r1     = basis1.primitive(j).center();

            auto gamma_ij = alpha0 + alpha1;
            auto rho_ij   = alpha0 * alpha1 / gamma_ij;
            auto dr2      = distance_squared(r0, r1);
            K[i][j]       = std::exp(-rho_ij * dr2);
        }
    }
    return K;
}

inline auto coarse_k_ij(const simde::type::ao_basis_set& basis0,
                        const simde::type::ao_basis_set& basis1) {
    auto K                   = k_ij(basis0, basis1);
    auto bs0_libint          = make_libint_basis_set(basis0);
    auto bs1_libint          = make_libint_basis_set(basis1);
    std::size_t prim0_offset = 0;
    for(const auto& shell0 : bs0_libint) {
        auto nprims0 = shell0.nprim();
        for(std::size_t prim0 = 0; prim0 < nprims0; ++prim0) {
            auto i                   = prim0_offset + prim0;
            auto coefi               = std::fabs(shell0.contr[0].coeff[prim0]);
            std::size_t prim1_offset = 0;
            for(const auto& shell1 : bs1_libint) {
                auto nprims1 = shell1.nprim();
                for(std::size_t prim1 = 0; prim1 < nprims1; ++prim1) {
                    auto j     = prim1_offset + prim1;
                    auto coefj = std::fabs(shell1.contr[0].coeff[prim1]);
                    K[i][j] *= coefi * coefj;
                }
                prim1_offset += nprims1;
            }
        }
        prim0_offset += nprims0;
    }
    return K;
}
} // namespace integrals::libint::detail_
