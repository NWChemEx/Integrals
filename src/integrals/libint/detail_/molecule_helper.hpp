#pragma once
#include <array>
#include <chemist/molecule/molecule.hpp>

namespace integrals::detail_ {

using LI_molecule =
  std::vector<std::pair<double, std::array<double, 3ul>>,
              std::allocator<std::pair<double, std::array<double, 3ul>>>>;

/** @brief Converts a chemist::Molecule object to a LibInt2 object
 *
 *  The LibInt2 molecule structure is store as a vector of pairs.
 *  Each pair consists of a double for the charge, and a double array
 *  of length 3 for the nuclear coordinates.
 *
 *  The chemist Molecule is a class holding a vector of atoms. The class
 *  provides an iterator for accessing the atoms.
 *
 *  @param[in] mol  The chemist Molecule object to be converted
 *  @returns        The molecule as a LibInt object
 */
inline auto make_libint_molecule(const chemist::Molecule& NWX_mol) {
    LI_molecule LI_mol{};
    for(auto NWX_atom = NWX_mol.begin(); NWX_atom != NWX_mol.end();
        NWX_atom++) {
        auto charge = double(NWX_atom->Z());
        auto coords = NWX_atom->coords();
        auto LI_atom =
          std::pair<double, std::array<double, 3ul>>(charge, coords);
        LI_mol.push_back(LI_atom);
    }
    return LI_mol;
}

} // namespace integrals::detail_
