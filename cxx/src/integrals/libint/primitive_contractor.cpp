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

#include "../utils/primitive_index_helpers.hpp"
#include "detail_/make_libint_basis_set.hpp"
#include "detail_/primitive_pair_estimators.hpp"
#include "libint.hpp"
#include <cmath>
#include <integrals/property_types.hpp>
#include <libint2/shell.h>
#include <limits>
#include <unsupported/Eigen/CXX11/Tensor>

namespace integrals::libint {
namespace {

const auto desc = R"(
Primitive Contractor
====================

This module computes four-center electron repulsion integrals (ERIs) in the
contracted AO basis by explicitly contracting primitive integrals with
renormalized contraction coefficients.

The algorithm proceeds in three steps:

1. The "Raw Primitive ERI4" submodule is called with the original braket.
   It internally decontracts each basis set (one shell per primitive,
   coefficient = 1.0) and returns the raw primitive ERI tensor with libint's
   automatic shell normalization disabled. The tensor has dimensions
   [n_prim_aos, n_prim_aos, n_prim_aos, n_prim_aos], where n_prim_aos counts
   all AO components across all decontracted shells (e.g. a p-shell with 3
   primitives contributes 9 decontracted AO indices: 3 primitives x 3
   components).

2. The "Primitive Normalization" submodule is called on each of the four
   contracted basis sets. It returns a vector of renormalized contraction
   coefficients c[i], one per decontracted AO index, in the same ordering as
   the decontracted basis (primitive index varying fastest within each
   (shell, ao_component) block).

3. Matrices of pair estimates are computed for bra and ket basis sets. K and
   K' are respectively coarse estimates for bra and ket basis sets. Q and Q' are
   respectively fine estimates for bra and ket basis sets.

4. The contracted AO integrals are formed by summing
   c0[i]*c1[j]*c2[k]*c3[l]*prim_ERI[i,j,k,l] with screening:

   Coarse (from `coarse_k_ij`):
   - Skip bra pair (i,j) if K[i,j] < t
   - Skip ket pair (k,l) if K'[k,l] < t
   - Skip quartet if K[i,j] * K'[k,l] <= t (libint keeps only sums > ln(t))
   - Skip if |Q[i,j] * Q'[k,l] / sqrt(gamma_bra + gamma_ket)| < t.
)";

} // namespace

using my_pt   = simde::ERI4;
using norm_pt = integrals::property_types::Normalize<simde::type::ao_basis_set>;

MODULE_CTOR(PrimitiveContractor) {
    satisfies_property_type<my_pt>();
    description(desc);
    add_submodule<my_pt>("Raw Primitive ERI4");
    add_submodule<norm_pt>("Primitive Normalization");
    add_input<double>("Screening Threshold")
      .set_default(1e-16)
      .set_description(
        "Coarse screening uses the PrimitivePairEstimator matrices; fine "
        "screening uses libint ShellPair geometry factors. Matches libint "
        "ScreeningMethod::Original.");
}

MODULE_RUN(PrimitiveContractor) {
    const auto& [braket] = my_pt::unwrap_inputs(inputs);
    const double thresh  = inputs.at("Screening Threshold").value<double>();

    auto bra = braket.bra();
    auto ket = braket.ket();

    const auto& bs0 = bra.first.ao_basis_set();
    const auto& bs1 = bra.second.ao_basis_set();
    const auto& bs2 = ket.first.ao_basis_set();
    const auto& bs3 = ket.second.ao_basis_set();

    // Step 1: get raw primitive ERIs (decontracted, unnormalized).
    auto& raw_mod      = submods.at("Raw Primitive ERI4");
    const auto& prim_T = raw_mod.run_as<my_pt>(braket);

    // Step 2: get renormalized contraction coefficients for each basis set.
    auto& norm_mod = submods.at("Primitive Normalization");
    const auto& c0 = norm_mod.run_as<norm_pt>(bs0);
    const auto& c1 = norm_mod.run_as<norm_pt>(bs1);
    const auto& c2 = norm_mod.run_as<norm_pt>(bs2);
    const auto& c3 = norm_mod.run_as<norm_pt>(bs3);

    // Step 3: coarse screening matrices from PrimitivePairEstimator.
    const auto K_bra = detail_::coarse_k_ij(bs0, bs1);
    const auto k_ket = detail_::coarse_k_ij(bs2, bs3);

    auto map0 = utils::build_prim_ao_to_cgto_map(bs0);
    auto map1 = utils::build_prim_ao_to_cgto_map(bs1);
    auto map2 = utils::build_prim_ao_to_cgto_map(bs2);
    auto map3 = utils::build_prim_ao_to_cgto_map(bs3);

    auto pmap0 = utils::build_prim_ao_to_prim_shell_map(bs0);
    auto pmap1 = utils::build_prim_ao_to_prim_shell_map(bs1);
    auto pmap2 = utils::build_prim_ao_to_prim_shell_map(bs2);
    auto pmap3 = utils::build_prim_ao_to_prim_shell_map(bs3);

    // Number of AOs per basis set
    std::array naos{bs0.n_aos(), bs1.n_aos(), bs2.n_aos(), bs3.n_aos()};

    // Number of primitive shells per basis set
    std::array nprims{c0.size(), c1.size(), c2.size(), c3.size()};

    using namespace tensorwrapper;
    const auto prim_data = buffer::get_raw_data<double>(prim_T.buffer());

    auto gamma_bra = detail_::gamma_ij(bs0, bs1);
    auto gamma_ket = detail_::gamma_ij(bs2, bs3);
    auto Q_bra     = detail_::fine_k_ij(bs0, bs1);
    auto Q_ket     = detail_::fine_k_ij(bs2, bs3);
    std::vector<double> ao_data(naos[0] * naos[1] * naos[2] * naos[3], 0.0);

    auto ao_offset = [&](std::size_t i, std::size_t j, std::size_t k,
                         std::size_t l) {
        return i * naos[1] * naos[2] * naos[3] + j * naos[2] * naos[3] +
               k * naos[3] + l;
    };
    auto prim_offset = [&](std::size_t i, std::size_t j, std::size_t k,
                           std::size_t l) {
        return i * nprims[1] * nprims[2] * nprims[3] +
               j * nprims[2] * nprims[3] + k * nprims[3] + l;
    };

    for(std::size_t pshell_i = 0; pshell_i < nprims[0]; ++pshell_i) {
        const auto ci = c0[pshell_i];
        const auto mu = map0[pshell_i];
        const auto pi = pmap0[pshell_i];

        for(std::size_t pshell_j = 0; pshell_j < nprims[1]; ++pshell_j) {
            const auto nu = map1[pshell_j];
            const auto pj = pmap1[pshell_j];

            const auto cij      = ci * c1[pshell_j];
            const auto K_ij     = K_bra[pi][pj];
            const auto gamma_ij = gamma_bra[pi][pj];
            const auto Q_ij     = Q_bra[pi][pj];

            if(K_ij < thresh) continue;
            for(std::size_t pshell_k = 0; pshell_k < nprims[2]; ++pshell_k) {
                const auto lam    = map2[pshell_k];
                const auto pk     = pmap2[pshell_k];
                const double cijk = cij * c2[pshell_k];

                for(std::size_t pshell_l = 0; pshell_l < nprims[3];
                    ++pshell_l) {
                    const auto cl  = c3[pshell_l];
                    const auto sig = map3[pshell_l];
                    const auto pl  = pmap3[pshell_l];

                    const double K_kl     = k_ket[pk][pl];
                    const double Q_kl     = Q_ket[pk][pl];
                    const double gamma_kl = gamma_ket[pk][pl];

                    if(K_kl < thresh) continue;
                    if(K_ij * K_kl <= thresh) continue;

                    auto Q_ijkl     = Q_ij * Q_kl;
                    auto gamma_ijkl = gamma_ij + gamma_kl;
                    auto pfac       = std::abs(Q_ijkl / std::sqrt(gamma_ijkl));
                    if(pfac < thresh) continue;

                    ao_data[ao_offset(mu, nu, lam, sig)] +=
                      cijk * cl *
                      prim_data[prim_offset(pshell_i, pshell_j, pshell_k,
                                            pshell_l)];
                }
            }
        }
    }

    tensorwrapper::shape::Smooth shape({naos[0], naos[1], naos[2], naos[3]});
    tensorwrapper::buffer::Contiguous buf(std::move(ao_data), shape);
    simde::type::tensor t(shape, std::move(buf));

    auto result = results();
    return my_pt::wrap_results(result, t);
}

} // namespace integrals::libint
