#pragma once
#include "integrals/types.hpp"
#include <libchemist/basis_set/ao_basis_set.hpp>
#include <libint2.hpp>

namespace nwx_libint {

template<typename T>
using NWX_basis = libchemist::AOBasisSet<T>;
using LI_basis  = libint2::BasisSet;
using size      = integrals::type::size;

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
std::vector<LI_basis> make_basis_sets(
  const std::vector<NWX_basis<double>>& sets);
std::vector<LI_basis> make_basis_sets(
  const std::vector<NWX_basis<float>>& sets);

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

/** @brief Find the shells that contain the specified AOs.
 *
 *  Given a LibInt2 BasisSet @p basis_set and a lower @p lower and upper @p
 * upper AO index, returns a std::vector of shell indices that contain the AOs
 * between
 *  @p lower and @p upper.
 *
 *  @param[in] basis_set The LibInt2 BasisSet containing the AOs
 *  @param[in] lower The lower value of the AO range
 *  @param[in] upper The upper value of the AO range
 *  @returns The std::vector of the shell indices
 */
std::vector<size> aos2shells(LI_basis basis_set, size lower, size upper);

} // namespace nwx_libint
