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
Primitive Normalization
=======================

Computes the normalization factor :math:`1/\sqrt{(p|p)}` for each primitive
in the provided basis set, where :math:`(p|p)` is the diagonal overlap
integral of the primitive with itself.

The basis set is first decontracted so that each primitive is treated as a
separate shell. A libint2 overlap engine is then constructed with
``embed_normalization_into_coefficients = false`` so that the raw Gaussian
exponents and coefficients are used without any renormalization. The diagonal
overlap :math:`(p|p)` is extracted for each primitive shell and its
reciprocal square root is returned.
)";

} // namespace

using decontract_pt = integrals::property_types::DecontractBasisSet;
using pt = integrals::property_types::Normalize<simde::type::ao_basis_set>;

MODULE_CTOR(PrimitiveNormalization) {
    satisfies_property_type<pt>();
    description(desc);

    add_submodule<decontract_pt>("Decontract Basis Set")
      .set_description("Module used to decontract the basis set into "
                       "individual primitives");
}

MODULE_RUN(PrimitiveNormalization) {
    const auto& [basis] = pt::unwrap_inputs(inputs);

    // Decontract so each primitive is its own shell
    auto& decontract_mod   = submods.at("Decontract Basis Set");
    const auto& prim_basis = decontract_mod.run_as<decontract_pt>(basis);

    // Build a libint basis set without embedding normalization into
    // coefficients
    auto libint_bs = detail_::make_libint_basis_set(prim_basis, false);

    const auto n_shells = libint_bs.size();

    // Build an overlap engine over the decontracted basis
    if(!libint2::initialized()) libint2::initialize();
    const auto max_nprims = libint2::max_nprim(libint_bs);
    const auto max_l      = libint2::max_l(libint_bs);
    libint2::Engine engine(libint2::Operator::overlap, max_nprims, max_l);
    engine.set_max_nprim(max_nprims);
    engine.set(libint2::BraKet::xs_xs);

    const auto& buf = engine.results();

    std::vector<double> norms;
    norms.reserve(prim_basis.n_primitives());

    for(std::size_t i = 0; i < n_shells; ++i) {
        engine.compute(libint_bs[i], libint_bs[i]);
        const auto* ints = buf[0];

        // Each decontracted shell has size() AOs; extract the diagonal elements
        const auto n_aos = libint_bs[i].size();
        for(std::size_t mu = 0; mu < n_aos; ++mu) {
            // Diagonal element of the (n_aos x n_aos) overlap block
            const auto diag = ints[mu * n_aos + mu];
            norms.push_back(1.0 / std::sqrt(diag));
        }
    }

    auto result = results();
    return pt::wrap_results(result, norms);
}

} // namespace integrals::libint
