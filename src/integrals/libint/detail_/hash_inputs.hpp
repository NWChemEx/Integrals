#pragma once
#include "type_traits.hpp"
#include <functional>
#include <libint2.hpp>
#include <simde/types.hpp>
#include <string>
#include <vector>

template<>
struct std::hash<std::vector<libint2::BasisSet>> {
    std::size_t operator()(
      const std::vector<libint2::BasisSet>& b) const noexcept {
        /// TODO: Expand this
        return std::hash<std::size_t>{}(b.size());
    }
};

namespace integrals::detail_ {

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

template<typename OpType>
inline std::size_t hash_operator(const OpType& op) {
    /// TODO: Expand this
    constexpr auto libint_op = integrals::op_v<OpType>;
    return std::hash<int>{}(int(libint_op));
}

/** @brief Hashes the inputs to an integral module to produce a fxn_id
 *
 *  @tparam OpType The type of the operator
 *  @param[in] bases The basis sets as a vector of libint basis sets
 *  @param[in] op The operator associated with the integral
 *  @param[in] thresh The threshold ofr integral precision
 *  @returns A hash as a string
 */
template<typename OpType>
std::string hash_inputs(const std::vector<libint2::BasisSet>& bases,
                        const OpType& op, double thresh) {
    std::size_t hash = hash_operator(op);
    combine_hash(hash, bases, thresh);
    return std::to_string(hash);
}

} // namespace integrals::detail_