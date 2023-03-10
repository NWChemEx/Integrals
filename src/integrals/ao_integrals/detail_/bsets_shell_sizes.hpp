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
#include <simde/types.hpp>

namespace integrals::ao_integrals::detail_ {

/** @brief Collect the sizes of the shells of the basis sets.
 *
 *
 *
 *  @param[in] bases A vector of basis sets.
 *  @returns The shell sizes for the basis sets as a vector of vectors.
 */
inline auto bsets_shell_sizes(
  const std::vector<simde::type::ao_basis_set>& bases) {
    auto N = bases.size();
    std::vector<std::vector<std::size_t>> shell_sizes(N);
    for(auto i = 0; i < N; ++i) {
        for(const auto& shell : bases[i].shells()) {
            shell_sizes[i].push_back(shell.size());
        }
    }
    return shell_sizes;
}

} // namespace integrals::ao_integrals::detail_