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

using simde::type::aos;
using simde::type::braket;
using simde::type::t_e_type;
using simde::type::v_ee_type;
using simde::type::v_en_type;

template<typename BraKetType>
TEMPLATED_MODULE_CTOR(AOIntegral, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;
    satisfies_property_type<my_pt>();
    description("Computes integrals with Libint");
}

template<typename BraKetType>
TEMPLATED_MODULE_RUN(AOIntegral, BraKetType) {
    using my_pt = simde::EvaluateBraKet<BraKetType>;
    auto rv     = results();
    return my_pt::wrap_results(rv);
}

#define TWO_INDEX_AOI(op) AOIntegral<braket<aos, op, aos>>
#define TEMPLATE_2INDEX(op) template struct TWO_INDEX_AOI(op)
#define ADD_2INDEX(op, key) mm.add_module<TWO_INDEX_AOI(op)>(key)

TEMPLATE_2INDEX(t_e_type);
TEMPLATE_2INDEX(v_en_type);
TEMPLATE_2INDEX(v_ee_type);

void ao_integrals_set_defaults(pluginplay::ModuleManager& mm) {
    // Set any default associations
}

void load_ao_integrals(pluginplay::ModuleManager& mm) {
    ADD_2INDEX(t_e_type, "Kinetic");
    ADD_2INDEX(v_en_type, "Nuclear");
    ADD_2INDEX(v_ee_type, "ERI2");
    ao_integrals_set_defaults(mm);
}

#undef TWO_INDEX_AOI
#undef TEMPLATE_2INDEX
#undef ADD_2INDEX

} // namespace integrals::ao_integrals