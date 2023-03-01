/*
 * Copyright 2023 NWChemEx-Project
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
#include "libint.hpp"
#include "make_libint_factory.hpp"

namespace integrals::libint {

#define ADD_FACTORY(N, op, key) \
    mm.add_module<MakeLibintFactory<N, op>>(key " Factory")

void load_libint_modules(pluginplay::ModuleManager& mm) {
    ADD_FACTORY(2, simde::type::el_el_coulomb, "ERI2");
    ADD_FACTORY(3, simde::type::el_el_coulomb, "ERI3");
    ADD_FACTORY(4, simde::type::el_el_coulomb, "ERI4");
    ADD_FACTORY(2, simde::type::el_kinetic, "Kinetic");
    ADD_FACTORY(2, simde::type::el_nuc_coulomb, "Nuclear");
    ADD_FACTORY(2, simde::type::el_identity, "Overlap");
    ADD_FACTORY(2, simde::type::el_el_stg, "STG2");
    ADD_FACTORY(3, simde::type::el_el_stg, "STG3");
    ADD_FACTORY(4, simde::type::el_el_stg, "STG4");
    ADD_FACTORY(2, simde::type::el_el_yukawa, "Yukawa2");
    ADD_FACTORY(3, simde::type::el_el_yukawa, "Yukawa3");
    ADD_FACTORY(4, simde::type::el_el_yukawa, "Yukawa4");
    ADD_FACTORY(2, simde::type::el_el_f12_commutator, "F12 2C");
    ADD_FACTORY(3, simde::type::el_el_f12_commutator, "F12 3C");
    ADD_FACTORY(4, simde::type::el_el_f12_commutator, "F12 4C");
    ADD_FACTORY(2, simde::type::el_dipole, "EDipole");
    ADD_FACTORY(2, simde::type::el_quadrupole, "EQuadrupole");
    ADD_FACTORY(2, simde::type::el_octupole, "EOctupole");
    ADD_FACTORY(4, simde::type::el_el_delta, "DOI");
}

#undef ADD_FACTORY

} // namespace integrals::libint
