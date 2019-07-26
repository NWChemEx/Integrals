#pragma once
#include <libchemist/ao_basis_set.hpp>
#include <libint2.hpp>

namespace nwx_libint {
/** @brief Converts a libchemist::BasisSet object @p bs to a LibInt2 BasisSet
 *         object
 *
 *  @param[in] bs The libchemist BasisSet to be converted
 *  @returns The basis set as a LibInt2 BasisSet object
 */
libint2::BasisSet make_basis(const libchemist::AOBasisSet& bs);
} // namespace nwx_libint
