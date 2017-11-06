#pragma once
#include <libint2.hpp>
#include <LibChemist/SetOfAtoms.hpp>


namespace nwx_libint{
/** @brief Converts a Pulsar BasisSet object \p bs to a LibInt2 BasisSet object
 *
 *  @returns The basis set as a LibInt2 BasisSet object
 */
libint2::BasisSet make_basis(const LibChemist::BasisSet& bs);
}
