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

namespace integrals::ao_integrals {

// --- Declare Integral Module Types -------------------------------------------
template<std::size_t N, typename OperatorType>
DECLARE_MODULE(AOIntegral);

// --- Forward External Template Declarations ----------------------------------

#define EXTERN_TEMPLATE(N, op) extern template struct AOIntegral<N, op>

EXTERN_TEMPLATE(2, simde::type::el_kinetic);

#undef EXTERN_TEMPLATE

// --- Declare Module Load Functions -------------------------------------------

void load_ao_integrals(pluginplay::ModuleManager& mm);
void ao_integrals_set_defaults(pluginplay::ModuleManager& mm);

} // namespace integrals::ao_integrals
