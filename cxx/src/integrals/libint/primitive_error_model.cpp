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
#include "detail_/primitive_pair_estimators.hpp"
#include "libint.hpp"
#include <cmath>
#include <integrals/integrals.hpp>
#include <stdexcept>
#include <string>

namespace integrals::libint {
namespace {

const auto desc = R"(
Primitive Error Model
=====================

Uncertainty tensor for ERI4 primitive screening aligned with
`PrimitiveContractor`: same decontracted primitive-AO quadruple loop and the
same coarse / fine gates as libint `ScreeningMethod::Original` (coarse
`coarse_k_ij`, fine `fine_k_ij` with `gamma_ij`).

For each decontracted quartet that would be skipped by the contractor, this
module accumulates one of three per-quartet estimates into the corresponding
contracted AO element: fixed `Tolerance`, coarse pair product
`K_ij * K_kl`, or fine metric `|Q_ij Q_kl| / sqrt(gamma_ij + gamma_kl)`.
)";

/** @return True iff `PrimitiveContractor` would `continue` (skip) this quartet.
 */
inline bool primitive_quartet_skipped(double K_ij, double K_kl, double Q_ij,
                                      double Q_kl, double gamma_ij,
                                      double gamma_kl, double thresh) {
    if(K_ij < thresh) return true;
    if(K_kl < thresh) return true;
    if(K_ij * K_kl <= thresh) return true;
    const double gamma_ijkl = gamma_ij + gamma_kl;
    const double pfac       = std::abs(Q_ij * Q_kl / std::sqrt(gamma_ijkl));
    if(pfac < thresh) return true;
    return false;
}

enum class ErrorEstimateKind { Tolerance, Coarse, Fine };

inline ErrorEstimateKind parse_error_estimate(const std::string& s) {
    if(s == "Tolerance") return ErrorEstimateKind::Tolerance;
    if(s == "Coarse") return ErrorEstimateKind::Coarse;
    if(s == "Fine") return ErrorEstimateKind::Fine;
    throw std::invalid_argument(
      "Primitive Error Model: \"Error estimate\" must be \"Tolerance\", "
      "\"Coarse\", or \"Fine\"");
}

inline double skip_increment(ErrorEstimateKind kind, double thresh, double K_ij,
                             double K_kl, double Q_ij, double Q_kl,
                             double gamma_ij, double gamma_kl) {
    switch(kind) {
        case ErrorEstimateKind::Tolerance: return thresh;
        case ErrorEstimateKind::Coarse: return K_ij * K_kl;
        case ErrorEstimateKind::Fine:
            return std::abs(Q_ij * Q_kl / std::sqrt(gamma_ij + gamma_kl));
    }
    return 0.0;
}

} // namespace

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::Uncertainty<eri4_pt>;

MODULE_CTOR(PrimitiveErrorModel) {
    satisfies_property_type<pt>();
    description(desc);

    add_input<std::string>("Error estimate")
      .set_default("Tolerance")
      .set_description(
        "Per skipped primitive quartet: \"Tolerance\" adds the screening "
        "threshold; \"Coarse\" adds K_ij*K_kl; \"Fine\" adds the fine-screen "
        "metric |Q_ij Q_kl|/sqrt(gamma_ij+gamma_kl).");
}

MODULE_RUN(PrimitiveErrorModel) {
    const auto& [braket, tol] = pt::unwrap_inputs(inputs);
    const auto kind =
      parse_error_estimate(inputs.at("Error estimate").value<std::string>());

    auto bra = braket.bra();
    auto ket = braket.ket();

    const auto& bs0 = bra.first.ao_basis_set();
    const auto& bs1 = bra.second.ao_basis_set();
    const auto& bs2 = ket.first.ao_basis_set();
    const auto& bs3 = ket.second.ao_basis_set();

    const auto K_bra     = detail_::coarse_k_ij(bs0, bs1);
    const auto K_ket     = detail_::coarse_k_ij(bs2, bs3);
    const auto gamma_bra = detail_::gamma_ij(bs0, bs1);
    const auto gamma_ket = detail_::gamma_ij(bs2, bs3);
    const auto Q_bra     = detail_::fine_k_ij(bs0, bs1);
    const auto Q_ket     = detail_::fine_k_ij(bs2, bs3);

    auto map0 = utils::build_prim_ao_to_cgto_map(bs0);
    auto map1 = utils::build_prim_ao_to_cgto_map(bs1);
    auto map2 = utils::build_prim_ao_to_cgto_map(bs2);
    auto map3 = utils::build_prim_ao_to_cgto_map(bs3);

    auto pmap0 = utils::build_prim_ao_to_prim_shell_map(bs0);
    auto pmap1 = utils::build_prim_ao_to_prim_shell_map(bs1);
    auto pmap2 = utils::build_prim_ao_to_prim_shell_map(bs2);
    auto pmap3 = utils::build_prim_ao_to_prim_shell_map(bs3);

    std::array<std::size_t, 4> naos{bs0.n_aos(), bs1.n_aos(), bs2.n_aos(),
                                    bs3.n_aos()};
    std::array<std::size_t, 4> nprims{map0.size(), map1.size(), map2.size(),
                                      map3.size()};

    using float_type = double;
    tensorwrapper::shape::Smooth shape({naos[0], naos[1], naos[2], naos[3]});
    std::vector<float_type> raw_buffer(shape.size(), 0.0);

    auto ao_offset = [&](std::size_t i, std::size_t j, std::size_t k,
                         std::size_t l) {
        return i * naos[1] * naos[2] * naos[3] + j * naos[2] * naos[3] +
               k * naos[3] + l;
    };

    for(std::size_t pshell_i = 0; pshell_i < nprims[0]; ++pshell_i) {
        const auto mu = map0[pshell_i];
        const auto pi = pmap0[pshell_i];

        for(std::size_t pshell_j = 0; pshell_j < nprims[1]; ++pshell_j) {
            const auto nu = map1[pshell_j];
            const auto pj = pmap1[pshell_j];

            const double K_ij     = K_bra[pi][pj];
            const double gamma_ij = gamma_bra[pi][pj];
            const double Q_ij     = Q_bra[pi][pj];

            for(std::size_t pshell_k = 0; pshell_k < nprims[2]; ++pshell_k) {
                const auto lam = map2[pshell_k];
                const auto pk  = pmap2[pshell_k];

                for(std::size_t pshell_l = 0; pshell_l < nprims[3];
                    ++pshell_l) {
                    const auto sig = map3[pshell_l];
                    const auto pl  = pmap3[pshell_l];

                    const double K_kl     = K_ket[pk][pl];
                    const double Q_kl     = Q_ket[pk][pl];
                    const double gamma_kl = gamma_ket[pk][pl];

                    if(!primitive_quartet_skipped(K_ij, K_kl, Q_ij, Q_kl,
                                                  gamma_ij, gamma_kl, tol)) {
                        continue;
                    }

                    const double inc = skip_increment(
                      kind, tol, K_ij, K_kl, Q_ij, Q_kl, gamma_ij, gamma_kl);
                    raw_buffer[ao_offset(mu, nu, lam, sig)] += inc;
                }
            }
        }
    }

    tensorwrapper::buffer::Contiguous buffer(std::move(raw_buffer), shape);
    simde::type::tensor error(shape, std::move(buffer));
    auto result = results();
    return pt::wrap_results(result, error);
}

} // namespace integrals::libint
