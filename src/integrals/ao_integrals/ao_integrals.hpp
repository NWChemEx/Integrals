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
#include <simde/simde.hpp>

namespace integrals::ao_integrals {

// -----------------------------------------------------------------------------
// -- Declare Module Load Functions
// -----------------------------------------------------------------------------

void load_ao_integrals(pluginplay::ModuleManager& mm);
void ao_integrals_set_defaults(pluginplay::ModuleManager& mm);

// -----------------------------------------------------------------------------
// -- Declare Integral Module Types
// -----------------------------------------------------------------------------

template<std::size_t N, typename OperatorType, bool direct>
DECLARE_MODULE(AOIntegral);

template<std::size_t L, typename OperatorType>
DECLARE_MODULE(AOIntegralMultipole);

DECLARE_MODULE(AOIntegralDOI);

// -----------------------------------------------------------------------------
// -- Forward External Template Declarations
// -----------------------------------------------------------------------------

#define EXTERN_INT_AND_DIRECT(N, op)                \
    extern template class AOIntegral<N, op, false>; \
    extern template class AOIntegral<N, op, true>

EXTERN_INT_AND_DIRECT(2, simde::type::el_el_coulomb);
EXTERN_INT_AND_DIRECT(3, simde::type::el_el_coulomb);
EXTERN_INT_AND_DIRECT(4, simde::type::el_el_coulomb);
EXTERN_INT_AND_DIRECT(2, simde::type::el_kinetic);
EXTERN_INT_AND_DIRECT(2, simde::type::el_nuc_coulomb);
EXTERN_INT_AND_DIRECT(2, simde::type::el_identity);
EXTERN_INT_AND_DIRECT(2, simde::type::el_el_stg);
EXTERN_INT_AND_DIRECT(3, simde::type::el_el_stg);
EXTERN_INT_AND_DIRECT(4, simde::type::el_el_stg);
EXTERN_INT_AND_DIRECT(2, simde::type::el_el_yukawa);
EXTERN_INT_AND_DIRECT(3, simde::type::el_el_yukawa);
EXTERN_INT_AND_DIRECT(4, simde::type::el_el_yukawa);
EXTERN_INT_AND_DIRECT(2, simde::type::el_el_f12_commutator);
EXTERN_INT_AND_DIRECT(3, simde::type::el_el_f12_commutator);
EXTERN_INT_AND_DIRECT(4, simde::type::el_el_f12_commutator);
EXTERN_INT_AND_DIRECT(4, simde::type::el_el_delta);

#undef EXTERN_INT_AND_DIRECT

extern template class AOIntegralMultipole<0, simde::type::el_dipole>;
extern template class AOIntegralMultipole<1, simde::type::el_quadrupole>;
extern template class AOIntegralMultipole<2, simde::type::el_octupole>;

} // namespace integrals::ao_integrals