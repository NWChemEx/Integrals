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
#include <integrals/property_types.hpp>
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

3. The contracted AO integrals are formed by:

     AO[mu, nu, lam, sig] =
       sum_{i in prims(mu)} sum_{j in prims(nu)}
       sum_{k in prims(lam)} sum_{l in prims(sig)}
         c0[i] * c1[j] * c2[k] * c3[l] * prim_ERI[i, j, k, l]

   where prims(mu) is the set of decontracted AO indices that belong to the
   same contracted AO mu (same shell, same AO component, all primitives).
)";

// Build a map from decontracted AO index -> contracted AO index.
// For contracted shell s with n_prims primitives and n_aos AO components,
// the decontracted indices (prim p, ao_component a) all map to the contracted
// AO index first_ao[s] + a.
// NOTE: bs.size() returns the number of centers, not shells.
//       bs.n_shells() returns the total number of shells (flattened).
std::vector<std::size_t> build_prim_to_ao_map(
  const simde::type::ao_basis_set& bs) {
    std::vector<std::size_t> map;
    std::size_t contracted_ao = 0;
    for(std::size_t s = 0; s < bs.n_shells(); ++s) {
        const auto& shell  = bs.shell(s);
        const auto n_prims = shell.n_primitives();
        const auto n_aos   = shell.size();
        for(std::size_t p = 0; p < n_prims; ++p) {
            for(std::size_t a = 0; a < n_aos; ++a) {
                map.push_back(contracted_ao + a);
            }
        }
        contracted_ao += n_aos;
    }
    return map;
}

} // namespace

using my_pt   = simde::ERI4;
using norm_pt = integrals::property_types::Normalize<simde::type::ao_basis_set>;

MODULE_CTOR(PrimitiveContractor) {
    satisfies_property_type<my_pt>();
    description(desc);
    add_submodule<my_pt>("Raw Primitive ERI4");
    add_submodule<norm_pt>("Primitive Normalization");
}

MODULE_RUN(PrimitiveContractor) {
    const auto& [braket] = my_pt::unwrap_inputs(inputs);

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

    // Build primitive -> contracted AO index maps
    auto map0 = build_prim_to_ao_map(bs0);
    auto map1 = build_prim_to_ao_map(bs1);
    auto map2 = build_prim_to_ao_map(bs2);
    auto map3 = build_prim_to_ao_map(bs3);

    const Eigen::Index n0 = bs0.n_aos();
    const Eigen::Index n1 = bs1.n_aos();
    const Eigen::Index n2 = bs2.n_aos();
    const Eigen::Index n3 = bs3.n_aos();

    const Eigen::Index np0 = c0.size();
    const Eigen::Index np1 = c1.size();
    const Eigen::Index np2 = c2.size();
    const Eigen::Index np3 = c3.size();

    // Get raw data from the primitive tensor via Eigen TensorMap
    using namespace tensorwrapper;
    const auto pdata = buffer::get_raw_data<double>(prim_T.buffer());

    using prim_map_t =
      Eigen::TensorMap<Eigen::Tensor<const double, 4, Eigen::RowMajor>>;
    prim_map_t prim(pdata.data(), np0, np1, np2, np3);

    // Step 3: contract into AO basis
    std::vector<double> ao_data(n0 * n1 * n2 * n3, 0.0);
    using ao_map_t =
      Eigen::TensorMap<Eigen::Tensor<double, 4, Eigen::RowMajor>>;
    ao_map_t ao(ao_data.data(), n0, n1, n2, n3);

    for(Eigen::Index i = 0; i < np0; ++i) {
        const double ci       = c0[i];
        const Eigen::Index mu = map0[i];
        for(Eigen::Index j = 0; j < np1; ++j) {
            const double cij      = ci * c1[j];
            const Eigen::Index nu = map1[j];
            for(Eigen::Index k = 0; k < np2; ++k) {
                const double cijk      = cij * c2[k];
                const Eigen::Index lam = map2[k];
                for(Eigen::Index l = 0; l < np3; ++l) {
                    const Eigen::Index sig = map3[l];
                    ao(mu, nu, lam, sig) += cijk * c3[l] * prim(i, j, k, l);
                }
            }
        }
    }

    // Wrap result in a TensorWrapper tensor
    tensorwrapper::shape::Smooth shape(
      {static_cast<std::size_t>(n0), static_cast<std::size_t>(n1),
       static_cast<std::size_t>(n2), static_cast<std::size_t>(n3)});
    tensorwrapper::buffer::Contiguous buf(std::move(ao_data), shape);
    simde::type::tensor t(shape, std::move(buf));

    auto result = results();
    return my_pt::wrap_results(result, t);
}

} // namespace integrals::libint
