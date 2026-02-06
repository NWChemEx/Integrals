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

/** @file property_types.hpp
 *
 * Property types defined by the Integrals library are defined here. It's likely
 * that this file will be moved into a distinct subdirectory if more property
 * types are added in the future.
 */
#pragma once
#include <simde/types.hpp>
#include <simde/utils/convert.hpp>

/** @namespace integrals::property_types
 *
 *  @brief The namespace for property types defined by the Integrals library
 */
namespace integrals::property_types {

// PT used to estimate the contribution of primitive pairs
DECLARE_PROPERTY_TYPE(PrimitivePairEstimator);
PROPERTY_TYPE_INPUTS(PrimitivePairEstimator) {
    using ao_basis = const simde::type::ao_basis_set&;
    auto rv        = pluginplay::declare_input()
                .add_field<ao_basis>("Bra Basis Set")
                .add_field<ao_basis>("Ket Basis Set");
    rv["Bra Basis Set"].set_description(
      "The atomic orbital basis set for the bra");
    rv["Ket Basis Set"].set_description(
      "The atomic orbital basis set for the ket");
    return rv;
}

PROPERTY_TYPE_RESULTS(PrimitivePairEstimator) {
    using tensor = simde::type::tensor;
    auto rv      = pluginplay::declare_result().add_field<tensor>(
      "Primitive Pair Estimates");
    rv["Primitive Pair Estimates"].set_description(
      "A tensor containing the estimated values for each primitive pair "
      "integral");
    return rv;
}

template<typename BasePT>
DECLARE_TEMPLATED_PROPERTY_TYPE(Uncertainty, BasePT);

template<typename BasePT>
TEMPLATED_PROPERTY_TYPE_INPUTS(Uncertainty, BasePT) {
    auto rv  = BasePT::inputs();
    auto rv0 = rv.template add_field<double>("Tolerance");
    rv0["Tolerance"].set_description("The screening threshold");
    return rv0;
}

template<typename BasePT>
TEMPLATED_PROPERTY_TYPE_RESULTS(Uncertainty, BasePT) {
    return BasePT::results();
}

using DecontractBasisSet =
  simde::Convert<simde::type::ao_basis_set, simde::type::ao_basis_set>;

} // end namespace integrals::property_types
