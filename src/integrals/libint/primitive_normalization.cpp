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
#include <integrals/property_types.hpp>

namespace integrals::libint {
namespace {

const auto desc = R"(
Primitive Normalization
=======================

Returns the renormalized contraction coefficients for each primitive in the
provided basis set, matching libint's internal ``renorm()`` convention.

For each contracted shell, libint applies two normalization steps:

1. A per-primitive factor :math:`N_p = \sqrt{2^l (2\zeta_p)^{l+3/2} /
   (\sqrt{\pi^3} (2l-1)!!)}` is multiplied into each raw coefficient
   :math:`d_p`.

2. All scaled coefficients are then divided by :math:`\sqrt{\langle\phi|\phi\rangle}`,
   where :math:`\langle\phi|\phi\rangle` is the self-overlap of the contracted
   shell computed with the already-scaled coefficients, so that the contracted
   shell has unit norm.

This module returns the resulting values :math:`d_p N_p / \sqrt{\langle\phi|\phi\rangle}`
by constructing the basis set via libint with ``embed_normalization = true``
and reading the coefficients back from the resulting shell objects.

The output vector has one entry per (primitive, AO component) pair, in the
same order as the decontracted basis: for each contracted shell, the
:math:`n_{\rm prims} \times n_{\rm AOs}` entries are listed with the primitive
index varying fastest.
)";

} // namespace

using pt = integrals::property_types::Normalize<simde::type::ao_basis_set>;

MODULE_CTOR(PrimitiveNormalization) {
    satisfies_property_type<pt>();
    description(desc);
}

MODULE_RUN(PrimitiveNormalization) {
    const auto& [basis] = pt::unwrap_inputs(inputs);

    // Build a libint basis set with normalization embedded into coefficients.
    // This triggers libint's renorm(), which applies both the per-primitive
    // normalization factor and the contracted-shell unit-norm scaling.
    auto libint_bs = detail_::make_libint_basis_set(basis, true);

    std::vector<double> norms;

    for(std::size_t s = 0; s < libint_bs.size(); ++s) {
        const auto& shell     = libint_bs[s];
        const auto n_prims    = shell.nprim();
        const auto& nwx_shell = basis.shell(s);
        const auto n_aos      = nwx_shell.size();

        // For each primitive, emit one entry per AO component. All AO
        // components within a primitive share the same renormalized coefficient
        // (libint stores one coefficient per primitive per contraction).
        for(std::size_t p = 0; p < n_prims; ++p) {
            // shell.contr[0].coeff[p] holds the renormalized coefficient after
            // libint's renorm(): d_p * N_p / sqrt(contracted-shell norm)
            const auto c = std::abs(shell.contr[0].coeff[p]);
            for(std::size_t ao = 0; ao < n_aos; ++ao) { norms.push_back(c); }
        }
    }

    auto result = results();
    return pt::wrap_results(result, norms);
}

} // namespace integrals::libint
