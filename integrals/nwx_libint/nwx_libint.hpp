#pragma once
#include "integrals/types.hpp"
#include <libint2.hpp>

namespace nwx_libint {
/** @brief Converts a libchemist::AOBasisSet object @p bs to a LibInt2 BasisSet
 *         object
 *
 *  @param[in] bs The libchemist BasisSet to be converted
 *  @returns The basis set as a LibInt2 BasisSet object
 */

template<typename T> using NWX_basis = integrals::type::basis_set<T>;
using LI_basis = libint2::BasisSet;
using size = integrals::type::size;

LI_basis make_basis(const NWX_basis<double>& bs);
LI_basis make_basis(const NWX_basis<float>& bs);

std::vector<LI_basis> make_basis_sets(const std::vector<NWX_basis<double>>& sets);
std::vector<LI_basis> make_basis_sets(const std::vector<NWX_basis<float>>& sets);

size sets_max_nprims(const std::vector<LI_basis>& sets);
int sets_max_l(const std::vector<LI_basis>& sets);

} // namespace nwx_libint
