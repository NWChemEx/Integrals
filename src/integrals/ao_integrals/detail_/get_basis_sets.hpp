/*
 * Copyright 2024 NWChemEx-Project
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
#include "make_libint_basis_set.hpp"
#include <simde/types.hpp>

namespace integrals::ao_integrals::detail_ {

/** @brief Deterimine how many basis sets are in the Bra and Ket.
 *
 *  @tparam BraType The type of the Bra
 *  @tparam KetType The type of the Ket
 *  @param[in] bra The Bra
 *  @param[in] ket The Ket
 *  @returns The number of basis sets in the Bra and Ket
 */
template<typename BraType, typename KetType>
constexpr int get_n(const BraType& bra, const KetType& ket) {
    using simde::type::aos;
    using simde::type::aos_squared;
    constexpr auto bra_is_aos         = std::is_same_v<BraType, aos>;
    constexpr auto bra_is_aos_squared = std::is_same_v<BraType, aos_squared>;
    constexpr auto ket_is_aos         = std::is_same_v<KetType, aos>;
    constexpr auto ket_is_aos_squared = std::is_same_v<KetType, aos_squared>;
    if constexpr(bra_is_aos && ket_is_aos) {
        return 2;
    } else if constexpr(bra_is_aos && ket_is_aos_squared) {
        return 3;
    } else if constexpr(bra_is_aos_squared && ket_is_aos_squared) {
        return 4;
    }
}

/** @brief Unpack the basis sets from the Bra and Ket.
 *
 *  @tparam BraType The type of the Bra
 *  @tparam KetType The type of the Ket
 *  @param[in] bra The Bra
 *  @param[in] ket The Ket
 *  @returns A vector of basis sets in libint form
 */
template<typename BraType, typename KetType>
std::vector<libint2::BasisSet> get_basis_sets(const BraType& bra,
                                              const KetType& ket) {
    using simde::type::aos;
    using simde::type::aos_squared;

    std::vector<libint2::BasisSet> basis_sets;

    if constexpr(std::is_same_v<BraType, aos>) {
        basis_sets.push_back(make_libint_basis_set(bra.ao_basis_set()));
    } else if constexpr(std::is_same_v<BraType, aos_squared>) {
        basis_sets.push_back(make_libint_basis_set(bra.lhs().ao_basis_set()));
        basis_sets.push_back(make_libint_basis_set(bra.rhs().ao_basis_set()));
    }

    if constexpr(std::is_same_v<KetType, aos>) {
        basis_sets.push_back(make_libint_basis_set(ket.ao_basis_set()));
    } else if constexpr(std::is_same_v<KetType, aos_squared>) {
        basis_sets.push_back(make_libint_basis_set(ket.lhs().ao_basis_set()));
        basis_sets.push_back(make_libint_basis_set(ket.rhs().ao_basis_set()));
    }

    return basis_sets;
}

} // namespace integrals::ao_integrals::detail_
