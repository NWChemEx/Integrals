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
#include "detail_/primitive_quartet_screener.hpp"
#include "libint.hpp"
#include <array>
#include <integrals/property_types.hpp>
#include <vector>

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

    // Step 3: the screening decisions and the primitive AO index maps.
    const detail_::PrimitiveQuartetScreener screener(bs0, bs1, bs2, bs3,
                                                     thresh);
    const utils::PrimitiveAOMaps maps0(bs0);
    const utils::PrimitiveAOMaps maps1(bs1);
    const utils::PrimitiveAOMaps maps2(bs2);
    const utils::PrimitiveAOMaps maps3(bs3);

    // Number of AOs per basis set
    std::array naos{bs0.n_aos(), bs1.n_aos(), bs2.n_aos(), bs3.n_aos()};

    // Number of primitive AOs per basis set
    std::array nprims{c0.size(), c1.size(), c2.size(), c3.size()};

    using namespace tensorwrapper;
    const auto prim_data = buffer::get_raw_data<double>(prim_T.buffer());

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

    using detail_::ScreeningReason;
    for(std::size_t pshell_i = 0; pshell_i < nprims[0]; ++pshell_i) {
        const auto ci = c0[pshell_i];
        const auto mu = maps0.to_cgto[pshell_i];
        const auto pi = maps0.to_prim_shell[pshell_i];

        for(std::size_t pshell_j = 0; pshell_j < nprims[1]; ++pshell_j) {
            const auto nu  = maps1.to_cgto[pshell_j];
            const auto pj  = maps1.to_prim_shell[pshell_j];
            const auto cij = ci * c1[pshell_j];

            if(screener.bra_pair_screened(pi, pj)) continue;
            for(std::size_t pshell_k = 0; pshell_k < nprims[2]; ++pshell_k) {
                const auto lam    = maps2.to_cgto[pshell_k];
                const auto pk     = maps2.to_prim_shell[pshell_k];
                const double cijk = cij * c2[pshell_k];

                for(std::size_t pshell_l = 0; pshell_l < nprims[3];
                    ++pshell_l) {
                    const auto sig = maps3.to_cgto[pshell_l];
                    const auto pl  = maps3.to_prim_shell[pshell_l];

                    if(screener.screen(pi, pj, pk, pl) != ScreeningReason::kept)
                        continue;

                    ao_data[ao_offset(mu, nu, lam, sig)] +=
                      cijk * c3[pshell_l] *
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
