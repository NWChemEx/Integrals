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

// #include "ao_integrals/ao_integrals.hpp"
// #include "libint/libint.hpp"
// #include "shapes/shapes.hpp"
// #include "transforms/transforms.hpp"
#include <integrals/integrals_mm.hpp>

// using namespace simde::type;

namespace integrals {

void set_defaults(pluginplay::ModuleManager& mm) {
    // /// Submodule name
    // auto fac_sub = "AO Integral Factory";

    // /// Set Factory and Shape for non-screened integrals
    // std::vector<std::string> module_names{
    //   "ERI2", "ERI3", "ERI4",    "Kinetic", "Nuclear", "Overlap", "STG2",
    //   "STG3", "STG4", "Yukawa2", "Yukawa3", "Yukawa4", "DOI4"};
    // for(const auto& name : module_names) {
    //     mm.change_submod(name, "Tensor Shape", "OneTileShape");
    //     mm.change_submod("Direct " + name, "Tensor Shape", "OneTileShape");
    //     mm.change_submod(name, fac_sub, name + " Factory");
    //     mm.change_submod("Direct " + name, fac_sub, name + " Factory");
    // }
    // mm.change_submod("STG 4 Center dfdr Squared", "Tensor Shape",
    //                  "OneTileShape");
    // mm.change_submod("STG 4 Center dfdr Squared", fac_sub, "F12 4C Factory");

    // /// Set Factory and Shape for screened integrals
    // module_names = {"ERI3", "ERI4", "Kinetic", "Nuclear", "Overlap",
    //                 "STG3", "STG4", "Yukawa3", "Yukawa4"};
    // for(const auto& name : module_names) {
    //     auto name_cs = name + " CS";
    //     mm.change_submod(name_cs, "Tensor Shape", "OneTileShape");
    //     mm.change_submod("Direct " + name_cs, "Tensor Shape", "OneTileShape");
    //     mm.change_submod(name_cs, fac_sub, name + " Factory");
    //     mm.change_submod("Direct " + name_cs, fac_sub, name + " Factory");
    // }

    // /// Set Factory and Shape for multipoles
    // module_names = {"EDipole", "EQuadrupole", "EOctupole"};
    // for(const auto& name : module_names) {
    //     mm.change_submod(name, fac_sub, name + " Factory");
    // }

    // /// Set Factory for shell norms
    // /// TODO: Uncomment after module optimization
    // // mm.change_submod("Shell Norms Overlap", fac_sub, "Overlap Factory");
    // // mm.change_submod("Shell Norms Coulomb", fac_sub, "ERI4 Factory");
    // // mm.change_submod("Shell Norms STG", fac_sub, "STG4 Factory");
    // // mm.change_submod("Shell Norms Yukawa", fac_sub, "Yukawa4 Factory");

    // /// TODO: Need to be removed. See module declarations.
    // std::string sh_norm = "Shell Norms";
    // mm.change_submod("Kinetic CS", sh_norm, sh_norm + " Overlap");
    // mm.change_submod("Nuclear CS", sh_norm, sh_norm + " Overlap");
    // mm.change_submod("Overlap CS", sh_norm, sh_norm + " Overlap");
    // mm.change_submod("ERI3 CS", sh_norm, sh_norm + " Coulomb");
    // mm.change_submod("ERI4 CS", sh_norm, sh_norm + " Coulomb");
    // mm.change_submod("STG3 CS", sh_norm, sh_norm + " STG");
    // mm.change_submod("STG4 CS", sh_norm, sh_norm + " STG");
    // mm.change_submod("Yukawa3 CS", sh_norm, sh_norm + " Yukawa");
    // mm.change_submod("Yukawa4 CS", sh_norm, sh_norm + " Yukawa");
    // mm.change_submod("Direct Kinetic CS", sh_norm, sh_norm + " Overlap");
    // mm.change_submod("Direct Nuclear CS", sh_norm, sh_norm + " Overlap");
    // mm.change_submod("Direct Overlap CS", sh_norm, sh_norm + " Overlap");
    // mm.change_submod("Direct ERI3 CS", sh_norm, sh_norm + " Coulomb");
    // mm.change_submod("Direct ERI4 CS", sh_norm, sh_norm + " Coulomb");
    // mm.change_submod("Direct STG3 CS", sh_norm, sh_norm + " STG");
    // mm.change_submod("Direct STG4 CS", sh_norm, sh_norm + " STG");
    // mm.change_submod("Direct Yukawa3 CS", sh_norm, sh_norm + " Yukawa");
    // mm.change_submod("Direct Yukawa4 CS", sh_norm, sh_norm + " Yukawa");
}

void load_modules(pluginplay::ModuleManager& mm) {
    // ao_integrals::load_ao_integrals(mm);
    // libint::load_libint_modules(mm);
    // transforms::load_transformed_integrals(mm);
    // shapes::load_modules(mm);

    // ao_integrals::ao_integrals_set_defaults(mm);
    set_defaults(mm);
}

} // namespace integrals
