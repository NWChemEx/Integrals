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

#include "ao_integrals/ao_integrals.hpp"
#include "libint/libint.hpp"
#include "utils/utils.hpp"
#include <integrals/integrals_mm.hpp>

namespace integrals {

PLUGIN_CTOR(Integrals) {
    set_name("integrals");

    // The defaults refer to modules from more than one subdirectory, so all
    // modules are added before any defaults are set
    ao_integrals::load_modules(modules());
    libint::load_modules(modules());
    utils::load_modules(modules());

    libint::set_defaults(modules());
    ao_integrals::set_defaults(modules());
    utils::set_defaults(modules());
    change_submod("AO integral driver", "Kinetic", "Kinetic");
    change_submod("AO integral driver", "Electron-Nuclear attraction",
                  "Nuclear");
    change_submod("Four center J builder", "Four-center ERI", "ERI4");
    change_submod("Four center K builder", "Four-center ERI", "ERI4");
    change_submod("Density Fitting Integral", "Three-center ERI", "ERI3");
    change_submod("Coulomb Metric", "Two-center ERI", "ERI2");
}

} // namespace integrals
