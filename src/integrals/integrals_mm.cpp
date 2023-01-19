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

#include "integrals/integrals_mm.hpp"
#include "libint/cs_screened_integrals.hpp"
#include "libint/libint.hpp"
#include "libint/shellnorms.hpp"
#include "transforms/transforms.hpp"

#define ADD_INT_WITH_DIRECT(N, op, key_base)       \
    mm.add_module<Libint<N, op, false>>(key_base); \
    mm.add_module<Libint<N, op, true>>("Direct " key_base)

#define ADD_CS_INT_WITH_DIRECT(N, op, key_base)      \
    mm.add_module<CSLibint<N, op, false>>(key_base); \
    mm.add_module<CSLibint<N, op, true>>("Direct " key_base)

using namespace simde::type;

namespace integrals {

template<std::size_t N, typename OpType>
void register_transformed_integral(pluginplay::ModuleManager& mm,
                                   std::string key) {
    using module_t = StandardTransform<N, OpType>;
    auto new_key   = "Transformed " + key;
    mm.add_module<module_t>(new_key);
    mm.change_submod(new_key, "integral kernel", key);
}

void load_libint_integrals(pluginplay::ModuleManager& mm) {
    mm.add_module<LibintDOI<false>>("DOI");
    mm.add_module<LibintDOI<true>>("Direct DOI");
    mm.add_module<LibintMultipole<0, el_dipole>>("EDipole");
    mm.add_module<LibintMultipole<1, el_quadrupole>>("EQuadrupole");
    mm.add_module<LibintMultipole<2, el_octupole>>("EOctupole");
    ADD_INT_WITH_DIRECT(2, el_el_coulomb, "ERI2");
    ADD_INT_WITH_DIRECT(3, el_el_coulomb, "ERI3");
    ADD_INT_WITH_DIRECT(4, el_el_coulomb, "ERI4");
    ADD_INT_WITH_DIRECT(2, el_kinetic, "Kinetic");
    ADD_INT_WITH_DIRECT(2, el_nuc_coulomb, "Nuclear");
    ADD_INT_WITH_DIRECT(2, el_identity, "Overlap");
    ADD_INT_WITH_DIRECT(2, el_el_stg, "STG2");
    ADD_INT_WITH_DIRECT(3, el_el_stg, "STG3");
    ADD_INT_WITH_DIRECT(4, el_el_stg, "STG4");
    ADD_INT_WITH_DIRECT(2, el_el_yukawa, "Yukawa2");
    ADD_INT_WITH_DIRECT(3, el_el_yukawa, "Yukawa3");
    ADD_INT_WITH_DIRECT(4, el_el_yukawa, "Yukawa4");
    ADD_CS_INT_WITH_DIRECT(2, el_kinetic, "Kinetic CS");
    ADD_CS_INT_WITH_DIRECT(2, el_nuc_coulomb, "Nuclear CS");
    ADD_CS_INT_WITH_DIRECT(2, el_identity, "Overlap CS");
    ADD_CS_INT_WITH_DIRECT(3, el_el_coulomb, "ERI3 CS");
    ADD_CS_INT_WITH_DIRECT(4, el_el_coulomb, "ERI4 CS");
    ADD_CS_INT_WITH_DIRECT(3, el_el_stg, "STG3 CS");
    ADD_CS_INT_WITH_DIRECT(4, el_el_stg, "STG4 CS");
    ADD_CS_INT_WITH_DIRECT(3, el_el_yukawa, "Yukawa3 CS");
    ADD_CS_INT_WITH_DIRECT(4, el_el_yukawa, "Yukawa4 CS");

    mm.add_module<ShellNormOverlap>("Shell Norms Overlap");
    mm.add_module<ShellNormCoulomb>("Shell Norms Coulomb");
    mm.add_module<ShellNormSTG>("Shell Norms STG");
    mm.add_module<ShellNormYukawa>("Shell Norms Yukawa");

    mm.change_submod("Kinetic CS", "Shell Norms", "Shell Norms Overlap");
    mm.change_submod("Nuclear CS", "Shell Norms", "Shell Norms Overlap");
    mm.change_submod("Overlap CS", "Shell Norms", "Shell Norms Overlap");
    mm.change_submod("ERI3 CS", "Shell Norms", "Shell Norms Coulomb");
    mm.change_submod("ERI4 CS", "Shell Norms", "Shell Norms Coulomb");
    mm.change_submod("STG3 CS", "Shell Norms", "Shell Norms STG");
    mm.change_submod("STG4 CS", "Shell Norms", "Shell Norms STG");
    mm.change_submod("Yukawa3 CS", "Shell Norms", "Shell Norms Yukawa");
    mm.change_submod("Yukawa4 CS", "Shell Norms", "Shell Norms Yukawa");
    mm.change_submod("Direct Kinetic CS", "Shell Norms", "Shell Norms Overlap");
    mm.change_submod("Direct Nuclear CS", "Shell Norms", "Shell Norms Overlap");
    mm.change_submod("Direct Overlap CS", "Shell Norms", "Shell Norms Overlap");
    mm.change_submod("Direct ERI3 CS", "Shell Norms", "Shell Norms Coulomb");
    mm.change_submod("Direct ERI4 CS", "Shell Norms", "Shell Norms Coulomb");
    mm.change_submod("Direct STG3 CS", "Shell Norms", "Shell Norms STG");
    mm.change_submod("Direct STG4 CS", "Shell Norms", "Shell Norms STG");
    mm.change_submod("Direct Yukawa3 CS", "Shell Norms", "Shell Norms Yukawa");
    mm.change_submod("Direct Yukawa4 CS", "Shell Norms", "Shell Norms Yukawa");
}

void load_transformed_libint_integrals(pluginplay::ModuleManager& mm) {
    // register_transformed_integral<pt::edipole<T>>(mm, "EDipole");
    // register_transformed_integral<pt::equadrupole<T>>(mm, "EQuadrupole");
    // register_transformed_integral<pt::eoctopole<T>>(mm, "EOctopole");
    // register_transformed_integral<pt::eri2c<T>>(mm, "ERI2");
    register_transformed_integral<3, el_el_coulomb>(mm, "ERI3");
    register_transformed_integral<4, el_el_coulomb>(mm, "ERI4");
    // register_transformed_integral<pt::kinetic<T>>(mm, "Kinetic");
    // register_transformed_integral<pt::nuclear<T>>(mm, "Nuclear");
    // register_transformed_integral<pt::overlap<T>>(mm, "Overlap");
    register_transformed_integral<2, el_kinetic>(mm, "Kinetic");
    register_transformed_integral<2, el_nuc_coulomb>(mm, "Nuclear");
}

void load_f12_integrals(pluginplay::ModuleManager& mm) {
    mm.add_module<Libint<4, el_el_f12_commutator, false>>(
      "STG 4 Center dfdr Squared");
}

void load_transformed_f12_integrals(pluginplay::ModuleManager& mm) {
    register_transformed_integral<4, el_el_f12_commutator>(
      mm, "STG 4 Center dfdr Squared");
    register_transformed_integral<4, el_el_stg>(mm, "STG4");
    register_transformed_integral<4, el_el_yukawa>(mm, "Yukawa4");
}

void load_misc_transforms(pluginplay::ModuleManager& mm) {
    mm.add_module<StandardTransform<2, el_scf_k>>("Transformed K");
    mm.add_module<StandardTransform<2, fock>>("Transformed Fock");
}

void load_modules(pluginplay::ModuleManager& mm) {
    load_libint_integrals(mm);
    load_transformed_libint_integrals(mm);
    load_f12_integrals(mm);
    load_transformed_f12_integrals(mm);
    load_misc_transforms(mm);
}

} // namespace integrals

#undef ADD_INT_WITH_DIRECT
#undef ADD_CS_INT_WITH_DIRECT