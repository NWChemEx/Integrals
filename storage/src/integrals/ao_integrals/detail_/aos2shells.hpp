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

/** @brief Find the shells that contain the specified AOs.
 *
 *  Given an vector of shell sizes @p shell_sizes and a lower @p lower and
 *  upper @p AO index, returns a std::vector of shell indices that contain the
 *  AOs between @p lower and @p upper.
 *
 *  @param[in] shell_sizes The size of the shells of the basis sets
 *  @param[in] lo The lower value of the AO range
 *  @param[in] up The upper value of the AO range
 *  @returns An std::vector of the shell indices
 */
inline auto aos2shells(const std::vector<std::size_t>& shell_sizes,
                       std::size_t lo, std::size_t up) {
    std::vector<std::size_t> return_vec;
    for(auto ishell = 0, offset = 0; ishell < shell_sizes.size(); ++ishell) {
        if(offset >= up) break;
        if(offset >= lo) return_vec.push_back(ishell);
        offset += shell_sizes[ishell];
    }
    return return_vec;
}

} // namespace integrals::ao_integrals::detail_
