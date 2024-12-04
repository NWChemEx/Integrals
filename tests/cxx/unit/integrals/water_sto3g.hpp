/*
 * Copyright 2024 NWChemEx-Project
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
#include <simde/simde.hpp>

namespace test {

inline simde::type::molecule water_molecule() {
    using atom_t     = simde::type::atom;
    using molecule_t = simde::type::molecule;
    atom_t O{"O", 8ul, 0.0, 0.0, -0.143222342980786, 0.0};
    atom_t H1{"H", 1ul, 0.0, 1.638033502034240, 1.136556880358410, 0.0};
    atom_t H2{"H", 1ul, 0.0, -1.638033502034240, 1.136556880358410, 0.0};
    return molecule_t{O, H1, H2};
}

inline simde::type::ao_basis_set water_sto3g_basis_set() {
    using ao_basis_t     = simde::type::ao_basis_set;
    using atomic_basis_t = simde::type::atomic_basis_set;
    using cg_t           = simde::type::contracted_gaussian;
    using point_t        = simde::type::point;
    using doubles_t      = std::vector<double>;

    auto mol   = water_molecule();
    point_t r0 = mol[0].as_nucleus();
    point_t r1 = mol[1].as_nucleus();
    point_t r2 = mol[2].as_nucleus();

    doubles_t cs0{0.15432897, 0.53532814, 0.44463454};
    doubles_t es0{130.7093200, 23.8088610, 6.4436083};
    doubles_t cs1{-0.09996723, 0.39951283, 0.70011547};
    doubles_t es1{5.0331513, 1.1695961, 0.3803890};
    doubles_t cs2{0.15591627, 0.60768372, 0.39195739};
    doubles_t es2{5.0331513, 1.1695961, 0.3803890};
    cg_t cg0(cs0.begin(), cs0.end(), es0.begin(), es0.end(), r0);
    cg_t cg1(cs1.begin(), cs1.end(), es1.begin(), es1.end(), r0);
    cg_t cg2(cs2.begin(), cs2.end(), es2.begin(), es2.end(), r0);
    atomic_basis_t o("sto-3g", 8, r0);
    o.add_shell(chemist::ShellType::pure, 0, cg0);
    o.add_shell(chemist::ShellType::pure, 0, cg1);
    o.add_shell(chemist::ShellType::pure, 1, cg2);

    doubles_t cs3{0.1543289673, 0.5353281423, 0.4446345422};
    doubles_t es3{3.425250914, 0.6239137298, 0.1688554040};
    cg_t cg3(cs3.begin(), cs3.end(), es3.begin(), es3.end(), r1);
    cg_t cg4(cs3.begin(), cs3.end(), es3.begin(), es3.end(), r2);
    atomic_basis_t h0("sto-3g", 1, r1);
    atomic_basis_t h1("sto-3g", 1, r2);
    h0.add_shell(chemist::ShellType::pure, 0, cg3);
    h1.add_shell(chemist::ShellType::pure, 0, cg4);

    ao_basis_t bs;
    bs.add_center(o);
    bs.add_center(h0);
    bs.add_center(h1);
    return bs;
}

} // namespace test
