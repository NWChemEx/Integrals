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
#include <array>
#include <chemist/enums.hpp> /// For ShellType
#include <libint2.hpp>
#include <simde/types.hpp>
#include <vector>

namespace integrals::ao_integrals::detail_ {

/** @brief Unpacks the basis sets from the inputs
 *
 *  @param[in] inputs The module inputs containing the basis sets.
 *  @returns A vector of the basis sets.
 */
template<std::size_t N, typename ModuleInputs>
auto unpack_bases(const ModuleInputs& inputs) {
    using ao_space_t = simde::type::ao_space;
    using ao_basis_t = simde::type::ao_basis_set;
    std::vector<ao_basis_t> aos(N);
    if constexpr(N == 2) {
        aos[0] = inputs.at("bra").template value<ao_space_t>().basis_set();
        aos[1] = inputs.at("ket").template value<ao_space_t>().basis_set();
    } else if constexpr(N == 3) {
        aos[0] = inputs.at("bra").template value<ao_space_t>().basis_set();
        aos[1] = inputs.at("ket 1").template value<ao_space_t>().basis_set();
        aos[2] = inputs.at("ket 2").template value<ao_space_t>().basis_set();
    } else if constexpr(N == 4) {
        aos[0] = inputs.at("bra 1").template value<ao_space_t>().basis_set();
        aos[1] = inputs.at("bra 2").template value<ao_space_t>().basis_set();
        aos[2] = inputs.at("ket 1").template value<ao_space_t>().basis_set();
        aos[3] = inputs.at("ket 2").template value<ao_space_t>().basis_set();
    }
    return aos;
}

} // namespace integrals::ao_integrals::detail_
