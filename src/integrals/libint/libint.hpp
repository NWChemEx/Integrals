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

#define EXTERN_INT_AND_DIRECT(N, op)            \
    extern template class Libint<N, op, false>; \
    extern template class Libint<N, op, true>

namespace integrals {

template<std::size_t N, typename OperatorType, bool direct>
DECLARE_MODULE(Libint);
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

template<bool direct>
DECLARE_MODULE(LibintDOI);
extern template class LibintDOI<false>;
extern template class LibintDOI<true>;

template<std::size_t L, typename OperatorType>
DECLARE_MODULE(LibintMultipole);

extern template class LibintMultipole<0, simde::type::el_dipole>;
extern template class LibintMultipole<1, simde::type::el_quadrupole>;
extern template class LibintMultipole<2, simde::type::el_octupole>;

} // namespace integrals

#undef EXTERN_INT_AND_DIRECT