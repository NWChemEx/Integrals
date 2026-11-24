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
#include <simde/evaluate_braket/evaluate_braket.hpp>
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

DECLARE_PROPERTY_TYPE(PairScreener);
PROPERTY_TYPE_INPUTS(PairScreener) {
    using ao_basis = const simde::type::ao_basis_set&;
    auto rv        = pluginplay::declare_input()
                .add_field<ao_basis>("Bra Basis Set")
                .add_field<ao_basis>("Ket Basis Set")
                .add_field<double>("Tolerance");
    return rv;
}

PROPERTY_TYPE_RESULTS(PairScreener) {
    using index_vector = std::vector<std::size_t>;
    using pair_vector  = std::vector<index_vector>;
    auto rv            = pluginplay::declare_result().add_field<pair_vector>(
      "Primitive Pairs Passing Screening");
    return rv;
}

// PT used to report how screening treats each primitive quartet of an ERI4
DECLARE_PROPERTY_TYPE(PrimitiveQuartetScreening);
PROPERTY_TYPE_INPUTS(PrimitiveQuartetScreening) {
    auto rv  = simde::ERI4::inputs();
    auto rv0 = rv.add_field<double>("Tolerance");
    rv0["Tolerance"].set_description("The screening threshold");
    return rv0;
}

PROPERTY_TYPE_RESULTS(PrimitiveQuartetScreening) {
    using tensor = simde::type::tensor;
    auto rv      = pluginplay::declare_result()
                .add_field<tensor>("Primitive Quartet Contributions")
                .add_field<tensor>("Screening Reasons")
                .add_field<tensor>("Coarse Estimates")
                .add_field<tensor>("Fine Estimates");
    rv["Primitive Quartet Contributions"].set_description(
      "For each primitive AO quartet, what it contributes to the contracted "
      "ERI, whether or not it is screened out. Indexed by primitive AO.");
    rv["Screening Reasons"].set_description(
      "For each primitive quartet, 0 if it is kept, otherwise the first "
      "screening gate it fails: 1 bra pair, 2 ket pair, 3 coarse product, 4 "
      "fine estimate. Indexed by primitive.");
    rv["Coarse Estimates"].set_description(
      "For each primitive quartet, the coarse estimate K_ij * K_kl. Indexed "
      "by primitive.");
    rv["Fine Estimates"].set_description(
      "For each primitive quartet, the fine estimate "
      "|Q_ij Q_kl| / sqrt(gamma_ij + gamma_kl). Indexed by primitive.");
    return rv;
}

using DecontractBasisSet =
  simde::Convert<simde::type::ao_basis_set, simde::type::ao_basis_set>;

template<typename T>
DECLARE_TEMPLATED_PROPERTY_TYPE(Normalize, T);

template<typename T>
TEMPLATED_PROPERTY_TYPE_INPUTS(Normalize, T) {
    auto rv =
      pluginplay::declare_input().add_field<const T&>("Object to Normalize");
    rv["Object to Normalize"].set_description("The object to normalize");
    return rv;
}

template<typename T>
TEMPLATED_PROPERTY_TYPE_RESULTS(Normalize, T) {
    auto rv = pluginplay::declare_result().add_field<std::vector<double>>(
      "Normalization Factors");
    rv["Normalization Factors"].set_description(
      "A vector of normalization factors, one per primitive");
    return rv;
}

} // end namespace integrals::property_types
