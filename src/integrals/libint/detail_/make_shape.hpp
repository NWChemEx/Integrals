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

namespace integrals::detail_ {

/** @brief Given a vector of basis sets, compute the shape of the corresponding
 *         integral.
 *
 *  @param[in] bases A vector of LibInt2 BasisSets
 *  @param[in] leading_extent An extent value added to the front of the shape
 *  @returns A unique_ptr for the resulting shape
 */
inline auto make_shape(const std::vector<libint2::BasisSet>& bases,
                       std::size_t leading_extent = 0) {
    using shape_t   = typename simde::type::tensor::shape_type;
    using extents_t = typename shape_t::extents_type;

    /// Count up the extents from the shells in the basis sets
    extents_t extents{};
    if(leading_extent != 0) extents.push_back(leading_extent);
    for(auto& set : bases) {
        std::size_t set_extent = 0;
        for(auto& shell : set) { set_extent += shell.size(); }
        extents.push_back(set_extent);
    }

    /// Make the unique pointer to the shape
    return std::make_unique<shape_t>(extents);
}

} // namespace integrals::detail_