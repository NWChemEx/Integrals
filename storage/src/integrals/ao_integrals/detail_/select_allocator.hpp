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
#include "hash_inputs.hpp"

namespace integrals::ao_integrals::detail_ {

/** @brief Choose type of and return an allocator
 *
 *  @tparam direct Whether the integral if direct or not
 *  @tparam T The type the tensors field
 *  @tparam OpType The type of the operator
 *  @param[in] bases The basis sets as a vector
 *  @param[in] op The operator associated with the integral
 *  @param[in] cs_thresh The Cauchy-Schwarz screening threshold
 *  @returns An allocator
 */
template<bool direct, typename T, typename OpType>
auto select_allocator(const std::vector<simde::type::ao_basis_set>& bases,
                      const OpType& op, double cs_thresh = -1.0) {
    if constexpr(direct) {
        auto fxn_id = hash_inputs(bases, op, cs_thresh);
        return tensorwrapper::tensor::allocator::direct_ta_allocator<T>(fxn_id);
    } else {
        return tensorwrapper::tensor::default_allocator<T>();
    }
}

} // namespace integrals::ao_integrals::detail_
