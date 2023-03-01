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

namespace integrals::libint {

// -----------------------------------------------------------------------------
// -- Declare MakeLibintFactory Module Type
// -----------------------------------------------------------------------------

template<typename std::size_t N, typename OperatorType>
DECLARE_MODULE(MakeLibintFactory);

// -----------------------------------------------------------------------------
// -- Forward External Template Declarations
// -----------------------------------------------------------------------------

#define EXTERN_DECLARE(N, op) extern template class MakeLibintFactory<N, op>

EXTERN_DECLARE(2, simde::type::el_el_coulomb);
EXTERN_DECLARE(3, simde::type::el_el_coulomb);
EXTERN_DECLARE(4, simde::type::el_el_coulomb);
EXTERN_DECLARE(2, simde::type::el_kinetic);
EXTERN_DECLARE(2, simde::type::el_nuc_coulomb);
EXTERN_DECLARE(2, simde::type::el_identity);
EXTERN_DECLARE(2, simde::type::el_el_stg);
EXTERN_DECLARE(3, simde::type::el_el_stg);
EXTERN_DECLARE(4, simde::type::el_el_stg);
EXTERN_DECLARE(2, simde::type::el_el_yukawa);
EXTERN_DECLARE(3, simde::type::el_el_yukawa);
EXTERN_DECLARE(4, simde::type::el_el_yukawa);
EXTERN_DECLARE(2, simde::type::el_el_f12_commutator);
EXTERN_DECLARE(3, simde::type::el_el_f12_commutator);
EXTERN_DECLARE(4, simde::type::el_el_f12_commutator);
EXTERN_DECLARE(2, simde::type::el_dipole);
EXTERN_DECLARE(2, simde::type::el_quadrupole);
EXTERN_DECLARE(2, simde::type::el_octupole);
EXTERN_DECLARE(4, simde::type::el_el_delta);

#undef EXTERN_DECLARE

} // namespace integrals::libint