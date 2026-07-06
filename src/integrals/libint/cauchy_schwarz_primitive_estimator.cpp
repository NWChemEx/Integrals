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
#include <wtf/wtf.hpp>

namespace integrals::libint {
namespace {

const auto desc = R"(
CauchySchwarz Primitive Pair Estimator
=======================================

For each primitive shell pair (pi, pj) from bra_basis x ket_basis, returns:

    cspe[pi][pj] = |c_pi * c_pj| * Q_CS(pi, pj)

where
  - c_pi = d_pi * N_pi / sqrt(shell_norm_i) is the libint-renormalized
    contraction coefficient (embed_normalization=true), identical to the
    PrimitiveNormalization output used by PrimitiveContractor.
  - Q_CS(pi, pj) = sqrt(max_{a,b} (pi_a pj_b | pi_a pj_b)_raw)
    is the shell-level Cauchy-Schwarz factor: the square root of the largest
    diagonal self-pair ERI over all angular-momentum components a of pi and b of
    pj, computed with libint normalization disabled so no N-factor enters.

By the ERI Cauchy-Schwarz inequality applied to specific AO components
(m_i, m_j):

    |(pi_{m_i} pj_{m_j} | pk_{m_k} pl_{m_l})_raw|
    ≤ sqrt((pi_{m_i} pj_{m_j} | pi_{m_i} pj_{m_j})_raw)
      * sqrt((pk_{m_k} pl_{m_l} | pk_{m_k} pl_{m_l})_raw)
    ≤ Q_CS(pi, pj) * Q_CS(pk, pl)

so the product cspe[pi][pj] * cspe[pk][pl] is a rigorous upper bound on the
magnitude of any single primitive-AO ERI contribution to a contracted AO element.
)";

} // namespace

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::PrimitivePairEstimator;

MODULE_CTOR(CauchySchwarzPrimitiveEstimator) {
    satisfies_property_type<pt>();
    description(desc);
    add_submodule<eri4_pt>("Raw Primitive ERI4");
}

MODULE_RUN(CauchySchwarzPrimitiveEstimator) {
    const auto&& [bra_basis, ket_basis] = pt::unwrap_inputs(inputs);

    const auto n_bra_prims = bra_basis.n_primitives();
    const auto n_ket_prims = ket_basis.n_primitives();

    // Compute the raw (unnormalized) self-pair ERIs for the primitive shells.
    // Raw Primitive ERI4 decontracts internally and disables libint
    // normalization, so prim4[a,b,a',b'] = raw_ERI without any N-factor.
    simde::type::aos bra_aos(bra_basis);
    simde::type::aos ket_aos(ket_basis);
    simde::type::aos_squared bra_pair(bra_aos, ket_aos);
    simde::type::v_ee_type v_ee{};
    chemist::braket::BraKet mnls(bra_pair, v_ee, bra_pair);
    const auto& prim4 = submods.at("Raw Primitive ERI4").run_as<eri4_pt>(mnls);

    // Contracted coefficients with libint normalization embedded (default):
    //   coeff[p] = d_p * N_p / sqrt(contracted_shell_norm)
    // This matches the PrimitiveNormalization module used by
    // PrimitiveContractor.
    auto bra_libint = detail_::make_libint_basis_set(bra_basis);
    auto ket_libint = detail_::make_libint_basis_set(ket_basis);

    using tensorwrapper::buffer::make_contiguous;
    const auto& eris = make_contiguous(prim4.buffer());

    using float_type = double;
    std::vector<float_type> data(n_bra_prims * n_ket_prims, 0.0);
    tensorwrapper::shape::Smooth shape({n_bra_prims, n_ket_prims});
    tensorwrapper::buffer::Contiguous buffer(std::move(data), shape);

    using iter_type    = std::size_t;
    using index_vector = std::vector<iter_type>;
    using wtf::fp::float_cast;

    index_vector abs_prim{0, 0};

    // Accumulated absolute AO offsets into the prim4 tensor dimensions.
    // For each contracted shell s with n_prims primitives and n_aos AO
    // components, primitive p occupies AOs [shell_ao_offset + p*n_aos,
    // shell_ao_offset + (p+1)*n_aos).
    std::size_t bra_prim_offset = 0;
    std::size_t bra_ao_offset   = 0;

    for(std::size_t s0 = 0; s0 < bra_basis.n_shells(); ++s0) {
        const auto& bra_shell = bra_libint.at(s0);
        assert(bra_shell.contr.size() == 1); // No general contraction support
        const auto& bra_coeff = bra_shell.contr[0].coeff;
        const auto n_prims_s0 = bra_coeff.size();
        const auto naos_s0    = bra_basis.shell(s0).size();

        for(std::size_t p0 = 0; p0 < n_prims_s0; ++p0) {
            const auto c_i    = std::fabs(bra_coeff[p0]);
            abs_prim[0]       = bra_prim_offset + p0;
            const auto off_pi = bra_ao_offset + p0 * naos_s0;

            std::size_t ket_prim_offset = 0;
            std::size_t ket_ao_offset   = 0;

            for(std::size_t s1 = 0; s1 < ket_basis.n_shells(); ++s1) {
                const auto& ket_shell = ket_libint.at(s1);
                assert(ket_shell.contr.size() == 1); // No general contractions
                const auto& ket_coeff = ket_shell.contr[0].coeff;
                const auto n_prims_s1 = ket_coeff.size();
                const auto naos_s1    = ket_basis.shell(s1).size();

                for(std::size_t p1 = 0; p1 < n_prims_s1; ++p1) {
                    const auto c_j    = std::fabs(ket_coeff[p1]);
                    abs_prim[1]       = ket_prim_offset + p1;
                    const auto off_pj = ket_ao_offset + p1 * naos_s1;

                    // Q_CS(pi,pj) = sqrt(max_{a,b} (pi_a pj_b | pi_a pj_b)_raw)
                    // This is the correct Cauchy-Schwarz factor: by the ERI
                    // positivity, (pi_a pj_b | pi_a pj_b) >= 0, and by C-S:
                    // |(pi_a pj_b | pk_c pl_d)| <= sqrt((pi_a pj_b | pi_a
                    // pj_b))
                    //                               * sqrt((pk_c pl_d | pk_c
                    //                               pl_d))
                    //                           <= Q_CS(pi,pj) * Q_CS(pk,pl).
                    double max_diag = 0.0;
                    for(iter_type a = 0; a < naos_s0; ++a) {
                        for(iter_type b = 0; b < naos_s1; ++b) {
                            index_vector idx4 = {off_pi + a, off_pj + b,
                                                 off_pi + a, off_pj + b};
                            const auto val =
                              float_cast<float_type>(eris.get_elem(idx4));
                            max_diag = std::max(max_diag, val);
                        }
                    }
                    buffer.set_elem(abs_prim, c_i * c_j * std::sqrt(max_diag));
                }
                ket_prim_offset += n_prims_s1;
                ket_ao_offset += n_prims_s1 * naos_s1;
            }
        }
        bra_prim_offset += n_prims_s0;
        bra_ao_offset += n_prims_s0 * naos_s0;
    }

    simde::type::tensor rv(shape, std::move(buffer));
    auto result = results();
    return pt::wrap_results(result, rv);
}

} // namespace integrals::libint
