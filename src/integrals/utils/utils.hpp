/*
 * Copyright 2025 NWChemEx-Project
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

/** @file utils.hpp
 *
 * Declaration of utility modules for the Integrals library, plus common
 * infrastructural functions.
 */
#pragma once
#include <pluginplay/pluginplay.hpp>
#include <simde/simde.hpp>

/** @namespace integrals::utils
 *
 *  @brief The namespace for general utility modules for integrals
 */
namespace integrals::utils {

DECLARE_MODULE(DecontractBasisSet);

inline void set_defaults(pluginplay::ModuleManager& mm) {}

inline void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<DecontractBasisSet>("Decontract Basis Set");
    set_defaults(mm);
}

} // end namespace integrals::utils
