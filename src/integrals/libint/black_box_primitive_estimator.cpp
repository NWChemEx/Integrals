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

const auto desc = R"(
Libint Black Box Primitive Pair Estimator
=========================================

This module computes the matrix :math:`K_{\mu\nu}` where :math:`\mu` and
:math:`\nu` index AOs in the decontracted basis (one AO per primitive per
angular momentum component), and:

.. math::

   K_{\mu\nu} = c_\mu c_\nu \exp\left(-\frac{\zeta_\mu \zeta_\nu}
   {\zeta_\mu + \zeta_\nu} |\mathbf{R}_\mu - \mathbf{R}_\nu|^2\right)

where :math:`c_\mu = |d_\mu| \cdot N_\mu` is the raw contraction coefficient
scaled by the primitive normalization factor :math:`N_\mu = 1/\sqrt{(p_\mu|p_\mu)}`
from the PrimitiveNormalization submodule.

N.B. The algorithm assumes that the bra and ket are different. If they are the
same, we can save time by using the fact that the matrix is symmetric.

)";

// Computes square of the distance between two NWX point-like objects
template<typename T>
auto distance_squared(T&& a, T&& b) {
    auto dx = a.x() - b.x();
    auto dy = a.y() - b.y();
    auto dz = a.z() - b.z();
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

using decontract_pt = integrals::property_types::DecontractBasisSet;
using normalize_pt =
  integrals::property_types::Normalize<simde::type::ao_basis_set>;
using pt = integrals::property_types::PrimitivePairEstimator;

MODULE_CTOR(BlackBoxPrimitiveEstimator) {
    satisfies_property_type<pt>();
    description(desc);
    // TODO: Add citation for Chemist paper

    add_submodule<normalize_pt>("Primitive Normalization")
      .set_description("Module used to compute per-AO primitive normalization "
                       "factors 1/sqrt((p|p))");
    add_submodule<decontract_pt>("Decontract Basis Set")
      .set_description(
        "Module used to decontract the basis set into "
        "individual primitives (provides exponents and centers)");
}

MODULE_RUN(BlackBoxPrimitiveEstimator) {
    const auto&& [bra, ket] = pt::unwrap_inputs(inputs);

    using iter_type  = std::size_t;
    using float_type = double; // TODO: Get from basis sets

    // Normalization factors: N_mu = 1/sqrt((p_mu|p_mu)), one per decontracted
    // AO
    auto& norm_mod        = submods.at("Primitive Normalization");
    const auto& bra_norms = norm_mod.run_as<normalize_pt>(bra);
    const auto& ket_norms = norm_mod.run_as<normalize_pt>(ket);

    // Decontracted bases provide exponents and centers (one shell per
    // primitive)
    auto& dec_mod         = submods.at("Decontract Basis Set");
    const auto& bra_prims = dec_mod.run_as<decontract_pt>(bra);
    const auto& ket_prims = dec_mod.run_as<decontract_pt>(ket);

    // K is indexed over AOs of the decontracted basis (one per angular
    // momentum component per primitive).
    const iter_type n_bra_aos = bra_norms.size();
    const iter_type n_ket_aos = ket_norms.size();

    tensorwrapper::shape::Smooth shape({n_bra_aos, n_ket_aos});
    std::vector<float_type> data(shape.size(), 0.0);
    tensorwrapper::buffer::Contiguous buffer(std::move(data), shape);

    iter_type abs_ao_b = 0; // Absolute AO index in bra decontracted basis

    // Iterate over contracted shells; within each, iterate over primitives.
    // The decontracted basis has one shell per primitive (in the same order),
    // so decontracted shell index sb maps to contracted shell s and primitive
    // pb.
    iter_type sb = 0; // index into bra_prims (decontracted)
    for(iter_type s = 0; s < bra.n_shells(); ++s) {
        const auto n_prims_b = bra.shell(s).n_primitives();
        const auto n_aos_b   = bra.shell(s).size();

        for(iter_type pb = 0; pb < n_prims_b; ++pb, ++sb) {
            const auto bra_shell = bra_prims.shell(sb);
            const auto zeta0     = bra_shell.primitive(0).exponent();
            const auto bra_ctr   = bra_shell.center();
            // Raw contraction coefficient from the original contracted basis
            const auto raw0 =
              std::fabs(bra.shell(s).primitive(pb).coefficient());

            iter_type abs_ao_k =
              0; // Absolute AO index in ket decontracted basis

            iter_type sk = 0; // index into ket_prims (decontracted)
            for(iter_type t = 0; t < ket.n_shells(); ++t) {
                const auto n_prims_k = ket.shell(t).n_primitives();
                const auto n_aos_k   = ket.shell(t).size();

                for(iter_type pk = 0; pk < n_prims_k; ++pk, ++sk) {
                    const auto ket_shell = ket_prims.shell(sk);
                    const auto zeta1     = ket_shell.primitive(0).exponent();
                    const auto ket_ctr   = ket_shell.center();
                    const auto raw1 =
                      std::fabs(ket.shell(t).primitive(pk).coefficient());

                    const auto dr2 = distance_squared(bra_ctr, ket_ctr);

                    for(iter_type ao_b = 0; ao_b < n_aos_b; ++ao_b) {
                        // c_mu = |d_mu| * N_mu  where N_mu =
                        // 1/sqrt((p_mu|p_mu))
                        const auto coeff0 = raw0 * bra_norms[abs_ao_b + ao_b];

                        for(iter_type ao_k = 0; ao_k < n_aos_k; ++ao_k) {
                            const auto coeff1 =
                              raw1 * ket_norms[abs_ao_k + ao_k];

                            // K bar: Eq. 11 in the SI of the Chemist paper
                            auto k01 =
                              compute_k(zeta0, zeta1, coeff0, coeff1, dr2);
                            std::vector<iter_type> idx{abs_ao_b + ao_b,
                                                       abs_ao_k + ao_k};
                            buffer.set_elem(idx, k01);
                        } // loop over ket AOs
                    } // loop over bra AOs

                    abs_ao_k += n_aos_k;
                } // loop over ket primitives
            } // loop over ket contracted shells

            abs_ao_b += n_aos_b;
        } // loop over bra primitives
    } // loop over bra contracted shells

    simde::type::tensor rv(shape, std::move(buffer));

    auto result = results();
    return pt::wrap_results(result, rv);
}

} // namespace integrals::libint
