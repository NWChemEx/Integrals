#pragma once
#include "type_traits.hpp"
#include <functional>
#include <libint2.hpp>
#include <simde/types.hpp>
#include <string>
#include <vector>

namespace integrals::detail_ {

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
    constexpr auto libint_op = integrals::op_v<OpType>;
    auto op_hash             = std::hash<int>{}(int(libint_op));
    auto thresh_hash         = std::hash<double>{}(thresh);
    auto bases_hash          = std::hash<std::size_t>{}(bases.size());

    auto hash = op_hash ^ (thresh_hash << 1);
    hash      = hash ^ (bases_hash << 1);
    return std::to_string(hash);
}

} // namespace integrals::detail_