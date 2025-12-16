/*
 * Copyright 2025 NWChemEx-Project
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

#include "utils.hpp"
#include <integrals/property_types.hpp>

namespace integrals::utils {

using pt          = integrals::property_types::decontract_basis_set;
using mol_basis_t = simde::type::ao_basis_set;
using abs_t       = simde::type::atomic_basis_set;
using cg_t        = simde::type::contracted_gaussian;
using doubles_t   = std::vector<double>;

MODULE_CTOR(DecontractBasisSet) {
    satisfies_property_type<pt>();
    description("Decontracts the input AO basis set");
}

MODULE_RUN(DecontractBasisSet) {
    const auto& [input_bs] = pt::unwrap_inputs(inputs);

    // Decontracted primitives all have a coefficient of 1.0.
    doubles_t coeffs{1.0};

    // Loop over the shells in each atomic basis set and decontract them.
    // For each primitive in a contracted Gaussian, make a new shell with
    // that only that primitive.
    mol_basis_t output_bs;
    for(const auto& abs : input_bs) {
        abs_t new_abs(abs.basis_set_name(), abs.atomic_number(),
                      abs.center().as_point());
        for(const auto& shell : abs) {
            for(const auto& prim : shell.contracted_gaussian()) {
                doubles_t exponent{prim.exponent()};
                cg_t cg(coeffs.begin(), coeffs.end(), exponent.begin(),
                        exponent.end(), abs.center().as_point());
                new_abs.add_shell(shell.pure(), shell.l(), cg);
            }
        }
        output_bs.add_center(new_abs);
    }

    auto result = results();
    return pt::wrap_results(result, output_bs);
}

} // end namespace integrals::utils
