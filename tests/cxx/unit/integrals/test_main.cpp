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

#include <libint2.hpp>
#include <vector>

namespace {

auto make_atoms() {
    using atom_type = libint2::Atom;
    std::vector<atom_type> atoms(4);

    atoms[0].atomic_number = 8;
    atoms[0].x             = 0.0;
    atoms[0].y             = -0.143222342980786;
    atoms[0].z             = 0.0;

    atoms[1].atomic_number = 8;
    atoms[1].x             = 3.0;
    atoms[1].y             = 3.0 - 0.143222342980786;
    atoms[1].z             = 3.0;

    atoms[2].atomic_number = 1;
    atoms[2].x             = 1.638033502034240;
    atoms[2].y             = 1.136556880358410;
    atoms[2].z             = 0.0;

    atoms[3].atomic_number = 1;
    atoms[3].x             = 4.638033502034240;
    atoms[3].y             = 4.136556880358410;
    atoms[3].z             = 3.0;

    return atoms;
}

/// Makes an atomic basis set with one s shell made of three primitives
auto make_atomic_basis() {
    using float_type         = double;
    using vector_type        = libint2::svector<float_type>;
    using contraction_type   = libint2::Shell::Contraction;
    using contraction_vector = libint2::svector<contraction_type>;
    using shell_type         = libint2::Shell;
    using atom_bases_type    = std::vector<shell_type>;

    const int l        = 0; // s shell
    const bool is_pure = false;
    vector_type coeffs{0.1543289673, 0.5353281423, 0.4446345422};
    contraction_type contraction{l, is_pure, coeffs};
    contraction_vector conts(1, contraction);

    vector_type exponents{3.425250914, 0.6239137298, 0.1688554040};

    auto atoms = make_atoms();
    atom_bases_type rv;
    for(const auto& atom : atoms) {
        std::array<float_type, 3> origin{atom.x, atom.y, atom.z};
        rv.push_back(shell_type(exponents, conts, origin));
    }
    return rv;
}

// auto make_basis_set() {
//     using atom_basis_type = std::vector<libint2::Shell>;
//     using basis_set_type  = libint2::BasisSet;

//     auto atoms = make_atoms();
//     auto abs   = make_atomic_basis();
//     std::vector<atom_basis_type> element_bases(4, abs);

//     return basis_set_type(atoms, element_bases);
// }

} // namespace

int main(int argc, char* argv[]) {
    libint2::initialize();
    libint2::Shell::do_enforce_unit_normalization(false);
    libint2::Engine engine(libint2::Operator::coulomb, 3, 0);
    const auto& buf_vec = engine.results();

    auto basis   = make_atomic_basis();
    auto& shell0 = basis[0];
    auto& shell1 = basis[1];
    auto& shell2 = basis[2];
    auto& shell3 = basis[3];
    std::cout << shell2.nprim() << std::endl;
    std::cout << shell2.alpha[0] << std::endl;
    std::cout << shell2.alpha[1] << std::endl;
    std::cout << shell2.alpha[2] << std::endl;
    std::cout << shell2.ncontr() << std::endl;
    std::cout << shell2.contr[0].coeff[0] << std::endl;
    std::cout << shell2.contr[0].coeff[1] << std::endl;
    std::cout << shell2.contr[0].coeff[2] << std::endl;

    // auto& shell3 = basis[3];

    engine.compute(shell0, shell1, shell2, shell3);
    if(buf_vec[0] == nullptr)
        throw std::runtime_error("Engine failed to compute integrals");

    std::cout << buf_vec[0][0] << std::endl;

    libint2::finalize();

    return 0;
}
