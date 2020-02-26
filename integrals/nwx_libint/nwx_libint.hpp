#pragma once
#include "integrals/types.hpp"
#include <libint2.hpp>

namespace nwx_libint {

template<typename T> using NWX_basis = integrals::type::basis_set<T>;
using LI_basis = libint2::BasisSet;
using size = integrals::type::size;

/** @brief Converts a libchemist::AOBasisSet object @p bs to a LibInt2 BasisSet
 *         object
 *
 *  @param[in] bs The libchemist BasisSet to be converted
 *  @returns The basis set as a LibInt2 BasisSet object
 */
LI_basis make_basis(const NWX_basis<double>& bs);
LI_basis make_basis(const NWX_basis<float>& bs);

/** @brief Converts a std::vector of libchemist::AOBasisSet objects @p sets to a
 *         std::vector of LibInt2 BasisSet objects
 *
 *  @param[in] sets The vector of libchemist basis sets to be converted
 *  @returns The vector of basis sets as a LibInt2 BasisSet objects
 */
std::vector<LI_basis> make_basis_sets(const std::vector<NWX_basis<double>>& sets);
std::vector<LI_basis> make_basis_sets(const std::vector<NWX_basis<float>>& sets);

/** @brief Given a std::vector of LibInt2 basis sets @p sets, find the maximum
 *         number of primitives in shells across of basis sets
 *
 *  @param[in] sets The vector of Libint2 basis sets to be checked
 *  @returns The largest number of primitives found in a shell
 */
size sets_max_nprims(const std::vector<LI_basis>& sets);

/** @brief Given a std::vector of LibInt2 basis sets @p sets, find the maximum
 *         angular momentum in shells across of basis sets
 *
 *  @param[in] sets The vector of Libint2 basis sets to be checked
 *  @returns The highest angular momentum found in a shell
 */
int sets_max_l(const std::vector<LI_basis>& sets);

} // namespace nwx_libint
