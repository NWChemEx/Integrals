/*
 * Copyright 2022 NWChemEx-Project
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
#include "libint_factory.hpp"
#include "libint_op.hpp"
#include "make_libint_factory.hpp"
#include <simde/integral_factory.hpp>

namespace integrals::libint {

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_CTOR(MakeLibintFactory, N, OperatorType) {
    using my_pt = simde::IntegralFactory<OperatorType>;

    satisfies_property_type<my_pt>();

    add_input<double>("Threshold")
      .set_default(1.0E-16)
      .set_description(
        "The target precision with which the integrals will be computed");
}

template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_RUN(MakeLibintFactory, N, OperatorType) {
    using my_pt               = simde::IntegralFactory<OperatorType>;
    using integral_factory    = simde::type::integral_factory;
    using libint_factory      = LibintFactory<N, OperatorType>;
    using libint_basis_vector = typename libint_factory::libint_basis_vector;

    /// Inputs
    const auto& [bases, op, deriv] = my_pt::unwrap_inputs(inputs);
    auto thresh                    = inputs.at("Threshold").value<double>();

    // Convert from NWX bases to libint
    libint_basis_vector libint_bases;
    for(const auto& basis_i : bases)
        libint_bases.push_back(detail_::make_libint_basis_set(basis_i));

    /// Make factory
    auto pfactory = std::make_unique<libint_factory>(std::move(libint_bases),
                                                     op, thresh, deriv);
    integral_factory fac(std::move(pfactory));
    auto rv = results();
    return my_pt::wrap_results(rv, fac);
}

// -----------------------------------------------------------------------------
// -- Template Declarations
// -----------------------------------------------------------------------------

#define TEMPLATE_DECLARE(N, op) template class MakeLibintFactory<N, op>

TEMPLATE_DECLARE(2, simde::type::el_el_coulomb);
TEMPLATE_DECLARE(3, simde::type::el_el_coulomb);
TEMPLATE_DECLARE(4, simde::type::el_el_coulomb);
TEMPLATE_DECLARE(2, simde::type::el_kinetic);
TEMPLATE_DECLARE(2, simde::type::el_nuc_coulomb);
TEMPLATE_DECLARE(2, simde::type::el_identity);
TEMPLATE_DECLARE(2, simde::type::el_el_stg);
TEMPLATE_DECLARE(3, simde::type::el_el_stg);
TEMPLATE_DECLARE(4, simde::type::el_el_stg);
TEMPLATE_DECLARE(2, simde::type::el_el_yukawa);
TEMPLATE_DECLARE(3, simde::type::el_el_yukawa);
TEMPLATE_DECLARE(4, simde::type::el_el_yukawa);
TEMPLATE_DECLARE(2, simde::type::el_el_f12_commutator);
TEMPLATE_DECLARE(3, simde::type::el_el_f12_commutator);
TEMPLATE_DECLARE(4, simde::type::el_el_f12_commutator);
TEMPLATE_DECLARE(2, simde::type::el_dipole);
TEMPLATE_DECLARE(2, simde::type::el_quadrupole);
TEMPLATE_DECLARE(2, simde::type::el_octupole);
TEMPLATE_DECLARE(4, simde::type::el_el_delta);

#undef TEMPLATE_DECLARE

} // namespace integrals::libint
