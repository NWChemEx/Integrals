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

#pragma once
#include "property_types/types.hpp"
#include <sde/property_type/property_type.hpp>

namespace property_types {

/**
 * @brief Property type for modules that build tensors with LibInt and
 * TiledArray
 *
 * @tparam ElementType The type of the element in the tensor. Defaults to
 * `double`.
 */
template<typename ElementType = double>
struct LibIntIntegral : public sde::PropertyType<LibIntIntegral<ElementType>> {
    /// Size type vector
    using size_vec = std::vector<type::size>;
    /// Pair vector
    using pair_vec = std::vector<std::pair<type::size, type::size>>;
    /// Generates the input fields required by this property type
    auto inputs_();
    /// Generates the result fields required by this property type
    auto results_();
}; // class OverlapIntegral

//------------------------Implementations---------------------------------------

template<typename ElementType>
auto LibIntIntegral<ElementType>::inputs_() {
    auto rv = sde::declare_input()
                .add_field<ElementType>("Threshold", 1.0E-16)
                .template add_field<size_vec>("Tile Size", size_vec{180})
                .template add_field<ElementType>("Screening Threshold", 0.0)
                .template add_field<pair_vec>("Atom Tile Groups", pair_vec{});
    rv["Threshold"].set_description("Convergence threshold of integrals");
    rv["Tile Size"].set_description(
      "Size threshold for tiling tensors by atom blocks");
    rv["Screening Threshold"].set_description(
      "Threshold for Cauchy-Schwarz screening");
    rv["Atom Tile Groups"].set_description(
      "Groups of Atoms to tile together. Overwrites Tile Size");
    return rv;
}

template<typename ElementType>
auto LibIntIntegral<ElementType>::results_() {
    auto rv = sde::declare_result();
    return rv;
}

extern template class LibIntIntegral<double>;
extern template class LibIntIntegral<float>;

} // namespace property_types
