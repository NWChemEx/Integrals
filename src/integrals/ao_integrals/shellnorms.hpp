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

namespace integrals::ao_integrals {

template<std::size_t NBodies, typename OperatorType>
DECLARE_MODULE(ShellNorms);

using ShellNormOverlap = ShellNorms<1, simde::type::el_identity>;
using ShellNormCoulomb = ShellNorms<2, simde::type::el_el_coulomb>;
using ShellNormSTG     = ShellNorms<2, simde::type::el_el_stg>;
using ShellNormYukawa  = ShellNorms<2, simde::type::el_el_yukawa>;

extern template class ShellNorms<1, simde::type::el_identity>;
extern template class ShellNorms<2, simde::type::el_el_coulomb>;
extern template class ShellNorms<2, simde::type::el_el_stg>;
extern template class ShellNorms<2, simde::type::el_el_yukawa>;

} // namespace integrals::ao_integrals
