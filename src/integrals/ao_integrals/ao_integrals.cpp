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

#include "ao_integrals.hpp"

namespace integrals::ao_integrals {

using chemist::Electron;
using chemist::qm_operator::Kinetic;

// --- Placeholder Elements ----------------------------------------------------
namespace placeholder {

#include <pluginplay/pluginplay.hpp>
template<typename OperatorType>
DECLARE_TEMPLATED_PROPERTY_TYPE(Placeholder, OperatorType);

template<typename OperatorType>
TEMPLATED_PROPERTY_TYPE_INPUTS(Placeholder, OperatorType) {
    auto rv = pluginplay::declare_input();
    return rv;
}

template<typename OperatorType>
TEMPLATED_PROPERTY_TYPE_RESULTS(Placeholder, OperatorType) {
    auto rv = pluginplay::declare_result();
    return rv;
}

} // namespace placeholder

// --- Define Module Constructor -----------------------------------------------
template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_CTOR(AOIntegral, N, OperatorType) {
    using my_pt = placeholder::Placeholder<OperatorType>;
    satisfies_property_type<my_pt>();
    description("Computes integrals with Libint");
}

// --- Define Module Run Function ----------------------------------------------
template<std::size_t N, typename OperatorType>
TEMPLATED_MODULE_RUN(AOIntegral, N, OperatorType) {
    using my_pt = placeholder::Placeholder<OperatorType>;
    auto rv = results();
    return my_pt::wrap_results(rv);
}

// --- Template Declarations ---------------------------------------------------
#define TEMPLATE_AOI(N, op) template struct AOIntegral<N, op>

TEMPLATE_AOI(2, Kinetic<Electron>);

#undef TEMPLATE_AOI

// --- Define Module Load Functions --------------------------------------------
#define ADD_AOI(N, op, key) mm.add_module<AOIntegral<N, op>>(key)

void ao_integrals_set_defaults(pluginplay::ModuleManager& mm) {
    // Set any default associations
}

void load_ao_integrals(pluginplay::ModuleManager& mm) {
    ADD_AOI(2, Kinetic<Electron>, "Kinetic");
    ao_integrals_set_defaults(mm);
}

#undef ADD_AOI

} // namespace integrals::ao_integrals