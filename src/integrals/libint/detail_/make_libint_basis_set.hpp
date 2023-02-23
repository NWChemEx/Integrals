/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once
#include <array>
#include <chemist/enums.hpp> /// For ShellType
#include <libint2.hpp>
#include <simde/types.hpp>
#include <vector>

namespace integrals::libint::detail_ {

/** @brief Converts an NWX basis set object to a LibInt2 basis set object.
 *
 *  @param[in] bs The NWX basis set to be converted
 *  @returns The basis set as a LibInt2 basis set
 */
inline auto make_libint_basis_set(const simde::type::ao_basis_set& bs) {
    /// Typedefs for everything
    using atom_t          = libint2::Atom;
    using shell_t         = libint2::Shell;
    using basis_t         = libint2::BasisSet;
    using cont_t          = libint2::Shell::Contraction;
    using svec_d_t        = libint2::svector<double>;
    using conts_t         = libint2::svector<cont_t>;
    using centers_t       = std::vector<atom_t>;
    using atom_bases_t    = std::vector<shell_t>;
    using element_bases_t = std::vector<atom_bases_t>;

    /// Inputs for BasisSet constructor
    centers_t centers{};
    element_bases_t element_bases{};

    /// Atom doesn't have a value ctor, so here's a stand in
    auto atom_ctor = [](int Z, double x, double y, double z) {
        atom_t atom{};
        atom.atomic_number = Z;
        atom.x             = x;
        atom.y             = y;
        atom.z             = z;
        return atom;
    };

    /// Origin for shell construction
    std::array<double, 3> origin = {0.0, 0.0, 0.0};

    /// Convert centers and their shells to libint equivalents.
    for(auto center_i = 0; center_i < bs.size(); ++center_i) {
        /// Add current center to atoms list
        auto& center = bs[center_i];
        centers.push_back(
          atom_ctor(center_i, center.x(), center.y(), center.z()));

        /// Gather shells for this center and add them to element_bases
        atom_bases_t atom_bases{};
        for(const auto&& shelli : center) {
            const auto nprims = shelli.n_unique_primitives();
            const auto prim0  = shelli.unique_primitive(0);
            const auto primN  = shelli.unique_primitive(nprims - 1);
            const bool pure   = shelli.pure() == chemist::ShellType::pure;
            const int l       = shelli.l();

            svec_d_t alphas(&prim0.exponent(), &primN.exponent() + 1);
            svec_d_t coefs(&prim0.coefficient(), &primN.coefficient() + 1);
            conts_t conts{cont_t{l, pure, coefs}};
            /// Use origin for position, because BasisSet moves shells to center
            atom_bases.push_back(shell_t(alphas, conts, origin));
        }
        element_bases.push_back(atom_bases);
    }

    /// Return the new basis set
    return basis_t(centers, element_bases);
}

}
