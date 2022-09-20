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

#pragma once
#include <pluginplay/module_base.hpp>
#include <simde/types.hpp>

namespace integrals {

template<std::size_t N, typename OperatorType, bool direct>
DECLARE_MODULE(CSLibint);
// Non-direct
extern template class CSLibint<3, simde::type::el_el_coulomb, false>;
extern template class CSLibint<4, simde::type::el_el_coulomb, false>;
extern template class CSLibint<3, simde::type::el_el_stg, false>;
extern template class CSLibint<4, simde::type::el_el_stg, false>;
extern template class CSLibint<3, simde::type::el_el_yukawa, false>;
extern template class CSLibint<4, simde::type::el_el_yukawa, false>;
// Direct
extern template class CSLibint<3, simde::type::el_el_coulomb, true>;
extern template class CSLibint<4, simde::type::el_el_coulomb, true>;
extern template class CSLibint<3, simde::type::el_el_stg, true>;
extern template class CSLibint<4, simde::type::el_el_stg, true>;
extern template class CSLibint<3, simde::type::el_el_yukawa, true>;
extern template class CSLibint<4, simde::type::el_el_yukawa, true>;

} // namespace integrals
