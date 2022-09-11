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
#include <chemist/chemist.hpp>
#include <libint2.hpp>

namespace integrals {
using NWX_basis    = chemist::AOBasisSet<double>;
using NWX_molecule = chemist::Molecule;
using LI_basis     = libint2::BasisSet;
/** @brief The property integrals we currently support
 */
enum class property {
    overlap, // Overlap integrals
    kinetic, // Kinetic energy integrals
    nuclear, // Nuclear attraction integrals
    eri      // Electron repulsion integrals
};
/** @brief A thin API to map NWX quantities to underlying integral libraries
 *
 *  Within NWChemEx there is a set up to deal with molecular structures,
 *  basis sets, shells, properties, and the like. Ultimately a library is used
 * to calculate the corresponding integrals. Currently we just use Libint2 but
 * in future it is possible that we might have multiple libraries to choose
 * from. In any case Libint2 has its own expectations for data types and the
 * like.
 *
 *  This class provides a simple interface that lets you create an integral
 *  factory instance, compute integrals, and clean up.
 */
class Factory {
private:
    static size_t instances;
    static NWX_molecule mol_null;
    std::unique_ptr<libint2::Engine> engine;
    LI_basis lbs1;
    LI_basis lbs2;
    LI_basis lbs3;
    LI_basis lbs4;
    libint2::Operator op;

public:
    Factory(const property prop, const NWX_basis& bs1, const NWX_basis& bs2,
            const NWX_basis& bs3, const NWX_basis& bs4,
            const NWX_molecule& mol);
    Factory(const property prop, const NWX_basis& bs1, const NWX_basis& bs2,
            const NWX_basis& bs3, const NWX_molecule& mol);
    Factory(const property prop, const NWX_basis& bs1, const NWX_basis& bs2,
            const NWX_molecule& mol);
    Factory(const property prop, const NWX_basis& bs1, const NWX_basis& bs2,
            const NWX_basis& bs3, const NWX_basis& bs4);
    Factory(const property prop, const NWX_basis& bs1, const NWX_basis& bs2,
            const NWX_basis& bs3);
    Factory(const property prop, const NWX_basis& bs1, const NWX_basis& bs2);
    const libint2::Engine::target_ptr_vec compute(const size_t s1,
                                                  const size_t s2);
    const libint2::Engine::target_ptr_vec compute(const size_t s1,
                                                  const size_t s2,
                                                  const size_t s3,
                                                  const size_t s4);
    ~Factory();
};
} // namespace integrals
