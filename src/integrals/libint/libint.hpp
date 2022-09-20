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

template<std::size_t N, typename OperatorType>
DECLARE_MODULE(Libint);

extern template class Libint<2, simde::type::el_el_coulomb>;
extern template class Libint<3, simde::type::el_el_coulomb>;
extern template class Libint<4, simde::type::el_el_coulomb>;
extern template class Libint<2, simde::type::el_kinetic>;
extern template class Libint<2, simde::type::el_nuc_coulomb>;
extern template class Libint<2, simde::type::el_identity>;
extern template class Libint<2, simde::type::el_el_stg>;
extern template class Libint<3, simde::type::el_el_stg>;
extern template class Libint<4, simde::type::el_el_stg>;
extern template class Libint<2, simde::type::el_el_yukawa>;
extern template class Libint<3, simde::type::el_el_yukawa>;
extern template class Libint<4, simde::type::el_el_yukawa>;
extern template class Libint<2, simde::type::el_el_f12_commutator>;
extern template class Libint<3, simde::type::el_el_f12_commutator>;
extern template class Libint<4, simde::type::el_el_f12_commutator>;

DECLARE_MODULE(LibintDOI);

template<std::size_t L, typename OperatorType>
DECLARE_MODULE(LibintMultipole);

extern template class LibintMultipole<0, simde::type::el_dipole>;
extern template class LibintMultipole<1, simde::type::el_quadrupole>;
extern template class LibintMultipole<2, simde::type::el_octupole>;

/// Direct Integrals
template<std::size_t N, typename OperatorType>
DECLARE_MODULE(LibintDirect);

extern template class LibintDirect<2, simde::type::el_el_coulomb>;
extern template class LibintDirect<3, simde::type::el_el_coulomb>;
extern template class LibintDirect<4, simde::type::el_el_coulomb>;
extern template class LibintDirect<2, simde::type::el_kinetic>;
extern template class LibintDirect<2, simde::type::el_nuc_coulomb>;
extern template class LibintDirect<2, simde::type::el_identity>;
extern template class LibintDirect<2, simde::type::el_el_stg>;
extern template class LibintDirect<3, simde::type::el_el_stg>;
extern template class LibintDirect<4, simde::type::el_el_stg>;
extern template class LibintDirect<2, simde::type::el_el_yukawa>;
extern template class LibintDirect<3, simde::type::el_el_yukawa>;
extern template class LibintDirect<4, simde::type::el_el_yukawa>;
extern template class LibintDirect<2, simde::type::el_el_f12_commutator>;
extern template class LibintDirect<3, simde::type::el_el_f12_commutator>;
extern template class LibintDirect<4, simde::type::el_el_f12_commutator>;

DECLARE_MODULE(LibintDirectDOI);

} // namespace integrals
