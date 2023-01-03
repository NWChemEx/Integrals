/*
 * Copyright 2023 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
#include <simde/simde.hpp>

namespace integrals::detail_ {

/** @brief Code factorization for building a tensor with AO basis sets.
 *
 *  The starting point for all transforms is the tensor in the AO basis set.
 *  To get the AO-based tensor we need to unpack the map of AO basis sets that
 *  the calling module has. This function wraps that process to avoid code
 *  duplication.
 *
 *  @tparam N The number of centers in the integral.
 *  @tparam OpType The type of the operator for the integral
 *
 *  @param[in] submod The module to use to compute the AO-based tensor.
 *  @param[in] all_aos The map from 0-based mode offset to the AO space to use
 *                     for it.
 *  @param[in] op      The operator in the integral.
 *
 *  @return The requested tensor.
 *
 *  @throw std::runtime_error if all_aos.size() != N. Strong throw guarantee.
 */
template<std::size_t N, typename OpType>
auto run_fundamental(pluginplay::Module& submod,
                     const simde::space_map_t<simde::type::ao_space>& all_aos,
                     const OpType& op) {
    using sub_pt = simde::AOTensorRepresentation<N, std::decay_t<OpType>>;

    if(all_aos.size() != N)
        throw std::runtime_error("One or modes doesn't have a basis set.");

    simde::type::tensor x;

    // There's at least two AO bases
    const auto& aos0 = all_aos.at(0).get();
    const auto& aos1 = all_aos.at(1).get();
    if constexpr(N == 2) {
        std::tie(x) = submod.run_as<sub_pt>(aos0, op, aos1);
    } else if constexpr(N == 3) {
        const auto& aos2 = all_aos.at(2).get();
        std::tie(x)      = submod.run_as<sub_pt>(aos0, op, aos1, aos2);
    } else if constexpr(N == 4) {
        const auto& aos2 = all_aos.at(2).get();
        const auto& aos3 = all_aos.at(3).get();
        std::tie(x)      = submod.run_as<sub_pt>(aos0, aos1, op, aos2, aos3);
    } else {
        static_assert(N == 2, "Unsupported number of centers");
    }
    return x;
}

} // namespace integrals::detail_
