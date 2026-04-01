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

#pragma once
#include <simde/simde.hpp>
#include <vector>

namespace integrals::utils {

/** @brief Build a map from a primitive AO index to a primitive shell index.
 *
 *  A contracted shell is a linear combination of primitive shells. Similarly,
 *  a contracted AO is a linear combination of primitive AOs. Given the offset
 *  of a primitive AO, this function returns the offset of the corresponding
 *  primitive shell.
 *
 *  @param[in] bs The basis set to build the map for.
 *  @return A vector whose length is the number of primitive AO indices in the
 *  basis set. The value at index i is the primitive shell index that the i-th
 *  primitive AO belongs to.
 */
inline std::vector<std::size_t> build_prim_ao_to_prim_shell_map(
  const simde::type::ao_basis_set& bs) {
    std::vector<std::size_t> map;
    std::size_t prim_offset = 0;
    for(std::size_t s = 0; s < bs.n_shells(); ++s) {
        const auto& shell  = bs.shell(s);
        const auto n_prims = shell.n_primitives();
        const auto n_aos   = shell.size();
        for(std::size_t p = 0; p < n_prims; ++p) {
            for(std::size_t a = 0; a < n_aos; ++a) {
                map.push_back(prim_offset + p);
            }
        }
        prim_offset += n_prims;
    }
    return map;
}

/** @brief Maps a primitive AO index to a contracted AO index.
 *
 *  @param[in] bs The basis set to build the map for.
 *  @return A vector whose length is the number of primitive AO indices in the
 *  basis set. The value at index i is the index of the contracted AO that
 *  the i-th primitive AO belongs to.
 */
inline std::vector<std::size_t> build_prim_ao_to_cgto_map(
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

} // namespace integrals::utils
