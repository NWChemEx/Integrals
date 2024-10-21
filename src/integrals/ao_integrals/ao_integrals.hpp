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
#include <pluginplay/pluginplay.hpp>
#include <simde/simde.hpp>

/** @namespace integrals::ao_integrals
 * 
 *  @brief The namespace for the modules that produce AO Integrals
 */
namespace integrals::ao_integrals {

/** @brief The Module for computing AO Integrals
 *  
 *  @tparam BraKetType The type of the BraKet input
 */
template<typename BraKetType>
DECLARE_MODULE(AOIntegral);

/** @brief Load the AO integral modules into a Module Manager
 *
 *  @param mm The Module Manager to load the modules into
 *
 *  @throw none No throw guarantee
 */
void load_ao_integrals(pluginplay::ModuleManager& mm);

/** @brief Set default module relationships
 *
 *  @param mm The Module Manager with modules whose defaults will be set
 *
 *  @throw none No throw guarantee
 */
void ao_integrals_set_defaults(pluginplay::ModuleManager& mm);

// Forward External Template Declarations
#define EXTERN_AOI extern template struct AOIntegral
#define BRAKET simde::type::braket
#define AOS simde::type::aos

EXTERN_AOI<BRAKET<AOS, simde::type::t_e_type, AOS>>;
EXTERN_AOI<BRAKET<AOS, simde::type::v_en_type, AOS>>;
EXTERN_AOI<BRAKET<AOS, simde::type::v_ee_type, AOS>>;

#undef EXTERN_AOI
#undef BRAKET
#undef AOS

} // namespace integrals::ao_integrals
