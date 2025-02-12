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

namespace integrals::ao_integrals {

DECLARE_MODULE(AOIntegralsDriver);
DECLARE_MODULE(JFourCenter);
DECLARE_MODULE(KFourCenter);

inline void set_defaults(pluginplay::ModuleManager& mm) {
    const auto ao_driver = "AO integral driver";
    mm.change_submod(ao_driver, "Coulomb matrix", "Four center J builder");
    mm.change_submod(ao_driver, "Exchange matrix", "Four center K builder");
}

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<AOIntegralsDriver>("AO integral driver");
    mm.add_module<JFourCenter>("Four center J builder");
    mm.add_module<KFourCenter>("Four center K builder");
    set_defaults(mm);
}

} // namespace scf::matrix_builder