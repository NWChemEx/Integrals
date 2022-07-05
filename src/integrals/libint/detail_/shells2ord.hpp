#pragma once
#include <libint2.hpp>

namespace integrals::detail_ {

/** @brief Find the ordinal indices spanned by the shell indices
 *
 *  @param[in] bases A vector of libint basis sets.
 *  @param[in] shell A vector of indices for shells in the basis sets.
 *  @returns An std::vector of the ordinal indices associated with the shells.
 */
inline auto shells2ord(const std::vector<libint2::BasisSet>& bases,
                       std::vector<std::size_t> shells) {
    /// Write me
    std::vector<std::size_t> ord_pos;
    return ord_pos;
}

} // namespace integrals::detail_