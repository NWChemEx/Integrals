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
#include <functional>
#include <simde/types.hpp>
#include <string>
#include <vector>

namespace integrals::ao_integrals::detail_ {

/** @brief Combines hashes in the same way as boost::hash_combine (apparently)
 *
 *  @note based on https://stackoverflow.com/a/2595226
 *
 *  @param[in,out] seed The hash we are modifying
 *  @param[in] addition The hash we are adding to the seed
 */
inline void hash_together(std::size_t& seed, const std::size_t& addition) {
    seed ^= addition + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/** @brief Default case for the corresponding variadic function
 *
 *  @note based on https://stackoverflow.com/a/2595226
 *
 *  @param[in,out] seed The hash we are modifying
 */
inline void combine_hash(std::size_t& seed) {}

/** @brief Hashes the input object and adds it to the seed hash
 *
 *  @note based on https://stackoverflow.com/a/2595226
 *
 *  @tparam T The type of the first object to hash
 *  @tparam Rest The types of the rest of the objects to hash
 *  @param[in,out] seed The hash we are modifying
 *  @param[in] v The object to hash and add to the seed
 *  @param[in] rest The other objects to hash and add
 */
template<typename T, typename... Rest>
inline void combine_hash(std::size_t& seed, const T& v, Rest... rest) {
    std::hash<T> hasher;
    hash_together(seed, hasher(v));
    combine_hash(seed, rest...);
}

/** @brief Produces a hash of the input operator
 *
 *  The libint_op and as_string should handle differentiating the operator
 *  types, but any additional information that can differentiate instances of
 *  an operator needs to be hashed as well.
 *
 *  This should eventually be replaced with std::hash specializations for the
 *  different operator types.
 *
 *  @tparam OpType The type of the operator we want to hash
 *  @param[in] op The operator that is being hashed
 *  @returns A hash of the operator
 */
template<typename OpType>
inline std::size_t hash_operator(const OpType& op) {
    /// NWX operator string
    auto hash = std::hash<std::string>{}(op.as_string());
    /// Additional information from op
    if constexpr(std::is_same_v<OpType, simde::type::el_nuc_coulomb>) {
        const auto& nuclei = op.template at<1>();
        for(const auto& ai : nuclei) {
            combine_hash(hash, ai.Z());
            for(const auto& coord : ai.coords()) combine_hash(hash, coord);
        }
    } else if constexpr(std::is_same_v<OpType, simde::type::el_el_stg> ||
                        std::is_same_v<OpType, simde::type::el_el_yukawa> ||
                        std::is_same_v<OpType,
                                       simde::type::el_el_f12_commutator>) {
        const auto& stg = op.template at<0>();
        combine_hash(hash, stg.coefficient);
        combine_hash(hash, stg.exponent);
    }
    return hash;
}

/** @brief Hashes the inputs to an integral module to produce a fxn_id
 *
 *  @tparam OpType The type of the operator
 *  @param[in] bases The basis sets as a vector of libint basis sets
 *  @param[in] op The operator associated with the integral
 *  @param[in] thresh The threshold for integral precision
 *  @param[in] cs_thresh The Cauchy-Schwarz screening threshold
 *  @returns A hash as a string
 */
template<typename OpType>
std::string hash_inputs(const std::vector<simde::type::ao_basis_set>& bases,
                        const OpType& op, double cs_thresh = -1.0) {
    std::size_t hash = hash_operator(op);
    combine_hash(hash, bases, cs_thresh);
    return std::to_string(hash);
}

} // namespace integrals::ao_integrals::detail_

/** @brief std::hash specializaton for a vector of LibInt basis sets */
template<>
struct std::hash<std::vector<simde::type::ao_basis_set>> {
    /** @brief Hashes the input vector
     *
     *  @param[in] b A vector of libint basis sets
     *  @returns A hash of the vector
     */
    std::size_t operator()(
      const std::vector<simde::type::ao_basis_set>& b) const noexcept {
        using integrals::ao_integrals::detail_::combine_hash;
        /// Start with the number of basis sets
        auto hash = std::hash<std::size_t>{}(b.size());
        /// Set through the sets and hash the info
        for(const auto& set : b) {
            for(const auto& p : set.unique_primitives()) {
                combine_hash(hash, p.coefficient(), p.exponent(), p.x(), p.y(),
                             p.z());
            }
        }
        return hash;
    }
};
