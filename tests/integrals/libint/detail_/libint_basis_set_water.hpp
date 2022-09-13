#pragma once
#include <libint2.hpp>

namespace testing {

/// Water STO-3G basis set in Libint Format
inline libint2::BasisSet water_basis_set() {
    using atom_t  = libint2::Atom;
    using shell_t = libint2::Shell;
    using basis_t = libint2::BasisSet;

    /// Constructor for atom_t
    auto atom_ctor = [](int Z, double x, double y, double z) {
        atom_t atom{};
        atom.atomic_number = Z;
        atom.x             = x;
        atom.y             = y;
        atom.z             = z;
        return atom;
    };

    /// Fill out atoms in water
    /// Z corresponds to Shells position in vector below,
    /// not actual Atomic number
    std::vector<atom_t> atoms{};
    atoms.push_back(atom_ctor(0, 0.0, -0.143222342980786, 0.0));
    atoms.push_back(atom_ctor(1, 1.638033502034240, 1.136556880358410, 0.0));
    atoms.push_back(atom_ctor(2, -1.638033502034240, 1.136556880358410, 0.0));

    /// Construct STO-3G shells, then group appropriately
    shell_t O_s1{{130.7093200, 23.8088610, 6.4436083},
                 {{0, true, {0.15432897, 0.53532814, 0.44463454}}},
                 {{0.0, -0.143222342980786, 0.0}}};
    shell_t O_s2{{5.0331513, 1.1695961, 0.3803890},
                 {{0, true, {-0.09996723, 0.39951283, 0.70011547}}},
                 {{0.0, -0.143222342980786, 0.0}}};
    shell_t O_p2{{5.0331513, 1.1695961, 0.3803890},
                 {{1, true, {0.15591627, 0.60768372, 0.39195739}}},
                 {{0.0, -0.143222342980786, 0.0}}};
    shell_t H1_s1{{3.425250914, 0.6239137298, 0.1688554040},
                  {{0, true, {0.1543289673, 0.5353281423, 0.4446345422}}},
                  {{1.638033502034240, 1.136556880358410, 0.0}}};
    shell_t H2_s1{{3.425250914, 0.6239137298, 0.1688554040},
                  {{0, true, {0.1543289673, 0.5353281423, 0.4446345422}}},
                  {{-1.638033502034240, 1.136556880358410, 0.0}}};
    std::vector<shell_t> O_shells{O_s1, O_s2, O_p2};
    std::vector<shell_t> H1_shells{H1_s1};
    std::vector<shell_t> H2_shells{H2_s1};
    std::vector<std::vector<shell_t>> shells{O_shells, H1_shells, H2_shells};

    /// Return the basis set
    return basis_t(atoms, shells);
}

} // namespace testing