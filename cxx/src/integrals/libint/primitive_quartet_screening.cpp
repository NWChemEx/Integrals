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

#include "detail_/primitive_quartet_screener.hpp"
#include "libint.hpp"
#include <array>
#include <integrals/property_types.hpp>
#include <vector>

namespace integrals::libint {
namespace {

const auto desc = R"(
Primitive Quartet Screening
===========================

This module reports, for every primitive quartet of a four-center electron
repulsion integral (ERI), whether libint's screening drops it and what it would
have contributed. It uses the same screening decisions as the "Primitive
Contractor ERI4" module, which matches libint's ``ScreeningMethod::Original``.

The results are:

- "Primitive Quartet Contributions". Element ``[i, j, k, l]`` is
  ``c_i c_j c_k c_l (ij|kl)``, where ``(ij|kl)`` comes from the "Raw Primitive
  ERI4" submodule and the ``c`` come from the "Primitive Normalization"
  submodule. Every quartet is included, whether it is screened out or not, so
  adding up the elements whose primitive quartet is screened out gives the
  screening error. The indices are primitive AO indices: for each contracted
  shell, its primitives in order, and for each primitive, the shell's AO
  components in order.
- "Screening Reasons". Element ``[p, q, r, s]`` is 0 if primitive quartet
  ``(pq|rs)`` is kept. Otherwise it is the first gate the quartet fails: 1 if
  the bra pair is screened out, 2 if the ket pair is screened out, 3 if the
  product of the coarse pair estimates is too small, and 4 if the fine
  estimate is too small. The indices are primitive indices: for each
  contracted shell, its primitives in order.
- "Coarse Estimates" and "Fine Estimates". The estimates the gates compare to
  the tolerance, indexed like "Screening Reasons".
)";

} // namespace

using pt      = integrals::property_types::PrimitiveQuartetScreening;
using eri4_pt = simde::ERI4;
using norm_pt = integrals::property_types::Normalize<simde::type::ao_basis_set>;

MODULE_CTOR(PrimitiveQuartetScreening) {
    satisfies_property_type<pt>();
    description(desc);
    add_submodule<eri4_pt>("Raw Primitive ERI4");
    add_submodule<norm_pt>("Primitive Normalization");
}

MODULE_RUN(PrimitiveQuartetScreening) {
    const auto& [braket, tol] = pt::unwrap_inputs(inputs);

    auto bra = braket.bra();
    auto ket = braket.ket();

    const auto& bs0 = bra.first.ao_basis_set();
    const auto& bs1 = bra.second.ao_basis_set();
    const auto& bs2 = ket.first.ao_basis_set();
    const auto& bs3 = ket.second.ao_basis_set();

    auto& raw_mod      = submods.at("Raw Primitive ERI4");
    const auto& prim_T = raw_mod.run_as<eri4_pt>(braket);

    auto& norm_mod = submods.at("Primitive Normalization");
    const auto& c0 = norm_mod.run_as<norm_pt>(bs0);
    const auto& c1 = norm_mod.run_as<norm_pt>(bs1);
    const auto& c2 = norm_mod.run_as<norm_pt>(bs2);
    const auto& c3 = norm_mod.run_as<norm_pt>(bs3);

    // Contributions, indexed by primitive AO
    using tensorwrapper::buffer::get_raw_data;
    const auto prim_data = get_raw_data<double>(prim_T.buffer());
    std::array n_prim_aos{c0.size(), c1.size(), c2.size(), c3.size()};
    std::vector<double> contributions(prim_data.size());
    std::size_t ijkl = 0;
    for(std::size_t i = 0; i < n_prim_aos[0]; ++i)
        for(std::size_t j = 0; j < n_prim_aos[1]; ++j)
            for(std::size_t k = 0; k < n_prim_aos[2]; ++k)
                for(std::size_t l = 0; l < n_prim_aos[3]; ++l, ++ijkl)
                    contributions[ijkl] =
                      c0[i] * c1[j] * c2[k] * c3[l] * prim_data[ijkl];

    // Screening decisions and estimates, indexed by primitive
    const detail_::PrimitiveQuartetScreener screener(bs0, bs1, bs2, bs3, tol);
    std::array n_prims{bs0.n_primitives(), bs1.n_primitives(),
                       bs2.n_primitives(), bs3.n_primitives()};
    const auto n_quartets = n_prims[0] * n_prims[1] * n_prims[2] * n_prims[3];
    std::vector<double> reasons(n_quartets);
    std::vector<double> coarse(n_quartets);
    std::vector<double> fine(n_quartets);
    std::size_t pqrs = 0;
    for(std::size_t p = 0; p < n_prims[0]; ++p)
        for(std::size_t q = 0; q < n_prims[1]; ++q)
            for(std::size_t r = 0; r < n_prims[2]; ++r)
                for(std::size_t s = 0; s < n_prims[3]; ++s, ++pqrs) {
                    const auto reason = screener.screen(p, q, r, s);
                    reasons[pqrs]     = static_cast<double>(reason);
                    coarse[pqrs]      = screener.coarse_estimate(p, q, r, s);
                    fine[pqrs]        = screener.fine_estimate(p, q, r, s);
                }

    using tensorwrapper::buffer::Contiguous;
    using tensorwrapper::shape::Smooth;
    auto make_tensor = [](std::vector<double> data, const auto& extents) {
        Smooth shape({extents[0], extents[1], extents[2], extents[3]});
        Contiguous buffer(std::move(data), shape);
        return simde::type::tensor(shape, std::move(buffer));
    };

    auto result = results();
    return pt::wrap_results(result,
                            make_tensor(std::move(contributions), n_prim_aos),
                            make_tensor(std::move(reasons), n_prims),
                            make_tensor(std::move(coarse), n_prims),
                            make_tensor(std::move(fine), n_prims));
}

} // namespace integrals::libint
