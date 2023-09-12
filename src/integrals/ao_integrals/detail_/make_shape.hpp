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
#include <libint2.hpp>
#include <simde/types.hpp>

namespace integrals::ao_integrals::detail_ {

/** @brief Given a vector of basis sets, compute the shape of the corresponding
 *         integral.
 *
 *  @param[in] bases A vector of basis sets
 *  @param[in] leading_extent An extent value added to the front of the shape
 *  @returns A unique_ptr for the resulting shape
 */
inline auto make_shape(const std::vector<simde::type::ao_basis_set>& bases,
                       std::size_t leading_extent = 0) {
    using shape_t   = typename simde::type::tensor::shape_type;
    using extents_t = typename shape_t::extents_type;

    /// Count up the extents from the shells in the basis sets
    extents_t extents{};
    if(leading_extent != 0) extents.push_back(leading_extent);
    for(const auto& set : bases) { extents.push_back(set.n_aos()); }

    /// Make the unique pointer to the shape
    return std::make_unique<shape_t>(extents);
}

} // namespace integrals::ao_integrals::detail_
