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

/** @namespace integrals::libint
 *
 *  @brief The namespace for the modules that produce AO Integrals
 */
namespace integrals::libint {

using simde::type::braket;

using simde::type::aos;
using simde::type::aos_squared;

using simde::type::op_base_type;
using simde::type::s_e_type;
using simde::type::t_e_type;
using simde::type::v_ee_type;
using simde::type::v_en_type;

/** @brief The Module for computing AO Integrals
 *
 *  @tparam BraKetType The type of the BraKet input
 */
template<typename BraKetType>
DECLARE_MODULE(Libint);

/** @brief Load the libint modules into a Module Manager
 *
 *  @param mm The Module Manager to load the modules into
 *
 *  @throw none No throw guarantee
 */
void load_libint(pluginplay::ModuleManager& mm);

/** @brief Set default module relationships
 *
 *  @param mm The Module Manager with modules whose defaults will be set
 *
 *  @throw none No throw guarantee
 */
void libint_set_defaults(pluginplay::ModuleManager& mm);

// Forward External Template Declarations
#define EXTERN_LIBINT extern template struct Libint

EXTERN_LIBINT<braket<aos, op_base_type, aos>>;
EXTERN_LIBINT<braket<aos, op_base_type, aos_squared>>;
EXTERN_LIBINT<braket<aos_squared, op_base_type, aos_squared>>;
EXTERN_LIBINT<braket<aos, s_e_type, aos>>;
EXTERN_LIBINT<braket<aos, t_e_type, aos>>;
EXTERN_LIBINT<braket<aos, v_en_type, aos>>;
EXTERN_LIBINT<braket<aos, v_ee_type, aos>>;
EXTERN_LIBINT<braket<aos, v_ee_type, aos_squared>>;
EXTERN_LIBINT<braket<aos_squared, v_ee_type, aos_squared>>;

#undef EXTERN_LIBINT

} // namespace integrals::libint
