#pragma once
#include "integrals/types.hpp"
#include <libint2.hpp>
#include <array>

namespace nwx_libint {
/** @brief Converts a libchemist::BasisSet object @p bs to a LibInt2 BasisSet
 *         object
 *
 *  @param[in] bs The libchemist BasisSet to be converted
 *  @returns The basis set as a LibInt2 BasisSet object
 */

template<typename T>
using coord_type = std::array<T, 3>;

template<typename T>
libint2::BasisSet make_basis(const integrals::type::basis_set<T>& bs);
} // namespace nwx_libint
