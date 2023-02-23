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
 *  Given an AOBasisSet @p basis_set and a lower @p lower and upper @p
 * AO index, returns a std::vector of shell indices that contain the AOs
 * between @p lower and @p upper.
 *
 *  @param[in] bs The AOBasisSet containing the AOs
 *  @param[in] lo The lower value of the AO range
 *  @param[in] up The upper value of the AO range
 *  @returns An std::vector of the shell indices
 */
inline auto aos2shells(const simde::ao_basis_set& bs, std::size_t lo,
                       std::size_t up) {
    std::vector<std::size_t> return_vec;
    auto shells = bs.shells();
    for(auto ishell = 0, offset = 0; ishell < shells.size(); ++ishell) {
        if(offset >= up) break;
        if(offset >= lo) return_vec.push_back(ishell);
        offset += shells[ishell].size();
    }
    return return_vec;
}

} // namespace integrals::ao_integrals::detail_
