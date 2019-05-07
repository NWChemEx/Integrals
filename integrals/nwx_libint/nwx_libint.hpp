#pragma once
#include <libint2.hpp>
#include <LibChemist/AOBasisSet.hpp>


namespace nwx_libint{
/** @brief Converts a LibChemist::BasisSet object @p bs to a LibInt2 BasisSet
 *         object
 *
 *  @param[in] bs The LibChemist BasisSet to be converted
 *  @returns The basis set as a LibInt2 BasisSet object
 */
libint2::BasisSet make_basis(const LibChemist::AOBasisSet& bs);
}
