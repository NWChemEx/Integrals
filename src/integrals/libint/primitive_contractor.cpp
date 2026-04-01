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

3. The "Primitive Pair Estimator" submodule supplies coarse screening only:
   bra matrix K[i,j] and ket matrix Q[k,l] (same layout as Black Box).

4. The contracted AO integrals are formed by summing
   c0[i]*c1[j]*c2[k]*c3[l]*prim_ERI[i,j,k,l] with screening:

   Coarse (from PrimitivePairEstimator, matches libint ln_scr / ShellPair
   pair pruning and the engine's ln_scr_bra + ln_scr_ket check):
   - Skip bra pair (i,j) if K[i,j] < t
   - Skip ket pair (k,l) if Q[k,l] < t
   - Skip quartet if K[i,j] * Q[k,l] <= t (libint keeps only sums > ln(t))

   Fine (libint engine prefactor check; uses ShellPair p.K and gammas, not
   the estimator matrix):
   - Skip if |Kgeom_bra * Kgeom_ket / sqrt(gamma_bra + gamma_ket) *
     c0*c1*c2*c3| < t, with Kgeom and gamma from libint2::ShellPair data
   (ShellPair built with a permissive ln_prec so all primitive pairs appear).
)";

} // namespace

using my_pt   = simde::ERI4;
using norm_pt = integrals::property_types::Normalize<simde::type::ao_basis_set>;
using est_pt  = integrals::property_types::PrimitivePairEstimator;

MODULE_CTOR(PrimitiveContractor) {
    satisfies_property_type<my_pt>();
    description(desc);
    add_submodule<my_pt>("Raw Primitive ERI4");
    add_submodule<norm_pt>("Primitive Normalization");
    add_submodule<est_pt>("Primitive Pair Estimator");
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
    auto& est_mod    = submods.at("Primitive Pair Estimator");
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

    const Eigen::Index n0 = bs0.n_aos();
    const Eigen::Index n1 = bs1.n_aos();
    const Eigen::Index n2 = bs2.n_aos();
    const Eigen::Index n3 = bs3.n_aos();

    const Eigen::Index np0 = c0.size();
    const Eigen::Index np1 = c1.size();
    const Eigen::Index np2 = c2.size();
    const Eigen::Index np3 = c3.size();

    using namespace tensorwrapper;
    const auto pdata = buffer::get_raw_data<double>(prim_T.buffer());

    using prim_map_t =
      Eigen::TensorMap<Eigen::Tensor<const double, 4, Eigen::RowMajor>>;
    prim_map_t prim(pdata.data(), np0, np1, np2, np3);

    const Eigen::Index nk0 = static_cast<Eigen::Index>(bs0.n_primitives());
    const Eigen::Index nk1 = static_cast<Eigen::Index>(bs1.n_primitives());
    const Eigen::Index nq2 = static_cast<Eigen::Index>(bs2.n_primitives());
    const Eigen::Index nq3 = static_cast<Eigen::Index>(bs3.n_primitives());

    // Fine screening: libint's p.K and gamma for every primitive pair
    // (permissive ln_prec so primpairs lists are complete).
    const double ln_all_pairs = std::numeric_limits<double>::lowest();
    auto bs0_libint           = detail_::make_libint_basis_set(bs0);
    auto bs1_libint           = detail_::make_libint_basis_set(bs1);
    auto bs2_libint           = detail_::make_libint_basis_set(bs2);
    auto bs3_libint           = detail_::make_libint_basis_set(bs3);

    std::vector<std::vector<double>> bra_geom_K(nk0,
                                                std::vector<double>(nk1, 0.0));
    std::vector<std::vector<double>> bra_gamma(nk0,
                                               std::vector<double>(nk1, 0.0));
    std::vector<std::vector<double>> ket_geom_K(nq2,
                                                std::vector<double>(nq3, 0.0));
    std::vector<std::vector<double>> ket_gamma(nq2,
                                               std::vector<double>(nq3, 0.0));

    {
        std::size_t abs_p0 = 0;
        for(std::size_t s0 = 0; s0 < bs0_libint.size(); ++s0) {
            const auto& sh0    = bs0_libint[s0];
            std::size_t abs_p1 = 0;
            for(std::size_t s1 = 0; s1 < bs1_libint.size(); ++s1) {
                const auto& sh1 = bs1_libint[s1];
                libint2::ShellPair sp;
                sp.init(sh0, sh1, ln_all_pairs,
                        libint2::ScreeningMethod::Original);
                for(const auto& pp : sp.primpairs) {
                    const std::size_t ai = abs_p0 + pp.p1;
                    const std::size_t aj = abs_p1 + pp.p2;
                    bra_geom_K[ai][aj]   = pp.K;
                    bra_gamma[ai][aj]    = sh0.alpha[pp.p1] + sh1.alpha[pp.p2];
                }
                abs_p1 += sh1.alpha.size();
            }
            abs_p0 += sh0.alpha.size();
        }
    }
    {
        std::size_t abs_p2 = 0;
        for(std::size_t s2 = 0; s2 < bs2_libint.size(); ++s2) {
            const auto& sh2    = bs2_libint[s2];
            std::size_t abs_p3 = 0;
            for(std::size_t s3 = 0; s3 < bs3_libint.size(); ++s3) {
                const auto& sh3 = bs3_libint[s3];
                libint2::ShellPair sp;
                sp.init(sh2, sh3, ln_all_pairs,
                        libint2::ScreeningMethod::Original);
                for(const auto& pp : sp.primpairs) {
                    const std::size_t ak = abs_p2 + pp.p1;
                    const std::size_t al = abs_p3 + pp.p2;
                    ket_geom_K[ak][al]   = pp.K;
                    ket_gamma[ak][al]    = sh2.alpha[pp.p1] + sh3.alpha[pp.p2];
                }
                abs_p3 += sh3.alpha.size();
            }
            abs_p2 += sh2.alpha.size();
        }
    }

    std::vector<double> ao_data(n0 * n1 * n2 * n3, 0.0);
    using ao_map_t =
      Eigen::TensorMap<Eigen::Tensor<double, 4, Eigen::RowMajor>>;
    ao_map_t ao(ao_data.data(), n0, n1, n2, n3);

    for(Eigen::Index i = 0; i < np0; ++i) {
        const double ci       = c0[i];
        const Eigen::Index mu = map0[i];
        const Eigen::Index pi = pmap0[i];
        for(Eigen::Index j = 0; j < np1; ++j) {
            const Eigen::Index pj = pmap1[j];
            const double K_ij     = K_bra[pi][pj];
            if(K_ij < thresh) continue;
            const double cij      = ci * c1[j];
            const Eigen::Index nu = map1[j];
            for(Eigen::Index k = 0; k < np2; ++k) {
                const double cijk      = cij * c2[k];
                const Eigen::Index lam = map2[k];
                const Eigen::Index pk  = pmap2[k];
                for(Eigen::Index l = 0; l < np3; ++l) {
                    const Eigen::Index pl = pmap3[l];
                    const double Q_kl     = k_ket[pk][pl];
                    if(Q_kl < thresh) continue;
                    if(K_ij * Q_kl <= thresh) continue;
                    const double Kgb  = bra_geom_K[pi][pj];
                    const double Kgk  = ket_geom_K[pk][pl];
                    const double gamb = bra_gamma[pi][pj];
                    const double gamk = ket_gamma[pk][pl];
                    const double pfac = std::abs(
                      Kgb * Kgk / std::sqrt(gamb + gamk) * cijk * c3[l]);
                    if(pfac < thresh) continue;
                    const Eigen::Index sig = map3[l];
                    ao(mu, nu, lam, sig) += cijk * c3[l] * prim(i, j, k, l);
                }
            }
        }
    }

    tensorwrapper::shape::Smooth shape(
      {static_cast<std::size_t>(n0), static_cast<std::size_t>(n1),
       static_cast<std::size_t>(n2), static_cast<std::size_t>(n3)});
    tensorwrapper::buffer::Contiguous buf(std::move(ao_data), shape);
    simde::type::tensor t(shape, std::move(buf));

    auto result = results();
    return my_pt::wrap_results(result, t);
}

} // namespace integrals::libint
