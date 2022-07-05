#pragma once
#include <libint2.hpp>

namespace integrals::detail_ {

/** @brief Find the shells that contain the specified AOs.
 *
 *  Given a LibInt2 BasisSet @p basis_set and a lower @p lower and upper @p
 * AO index, returns a std::vector of shell indices that contain the AOs
 * between @p lower and @p upper.
 *
 *  @param[in] bs The LibInt2 BasisSet containing the AOs
 *  @param[in] lo The lower value of the AO range
 *  @param[in] up The upper value of the AO range
 *  @returns An std::vector of the shell indices
 */
inline auto aos2shells(const libint2::BasisSet& bs, std::size_t lo,
                       std::size_t up) {
    std::vector<std::size_t> return_vec;
    for(auto ishell = 0, offset = 0; ishell < bs.size(); ++ishell) {
        if(offset >= up) break;
        if(offset >= lo) return_vec.push_back(ishell);
        offset += bs[ishell].size();
    }
    return return_vec;
}

} // namespace integrals::detail_