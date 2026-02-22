#pragma once
#include <simde/simde.hpp>

namespace integrals::testing {

// Inputs for H2 tests
inline simde::type::molecule h2_molecule() {
    using atom_t     = simde::type::atom;
    using molecule_t = simde::type::molecule;
    atom_t H0("H", 1ul, 1836.15, 0.0, 0.0, 0.0);
    atom_t H1("H", 1ul, 1836.15, 0.0, 0.0, 1.3984);
    return molecule_t{H0, H1};
}

// Inputs for Water tests
inline simde::type::molecule water_molecule() {
    using atom_t     = simde::type::atom;
    using molecule_t = simde::type::molecule;
    atom_t O{"O", 8ul, 0.0, 0.0, -0.143222342980786, 0.0};
    atom_t H1{"H", 1ul, 0.0, 1.638033502034240, 1.136556880358410, 0.0};
    atom_t H2{"H", 1ul, 0.0, -1.638033502034240, 1.136556880358410, 0.0};
    return molecule_t{O, H1, H2};
}

} // namespace integrals::testing