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
#include <integrals/integrals_mm.hpp>

namespace integrals {

/** @brief Set default module relationships
 *
 *  @param mm The Module Manager with modules whose defaults will be set
 *
 *  @throw none No throw guarantee
 */
void set_defaults(pluginplay::ModuleManager& mm) {
    mm.change_submod("AO integral driver", "Kinetic", "Kinetic");
    mm.change_submod("AO integral driver", "Electron-Nuclear attraction",
                     "Nuclear");
    mm.change_submod("Four center J builder", "Four-center ERI", "ERI4");
    mm.change_submod("Four center K builder", "Four-center ERI", "ERI4");
    mm.change_submod("Density Fitting Integral", "Three-center ERI", "ERI3");
    mm.change_submod("Coulomb Metric", "Two-center ERI", "ERI2");
}

void load_modules(pluginplay::ModuleManager& mm) {
    ao_integrals::load_modules(mm);
    libint::load_modules(mm);
    set_defaults(mm);
}

} // namespace integrals
