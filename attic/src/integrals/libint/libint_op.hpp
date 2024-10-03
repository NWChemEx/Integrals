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
#include <libint2.hpp>
#include <simde/types.hpp>

namespace integrals::libint {
namespace detail_ {

template<typename T>
struct LibintOp;

template<>
struct LibintOp<simde::type::el_el_delta> {
    static constexpr auto value = libint2::Operator::delta;
};

template<>
struct LibintOp<simde::type::el_el_coulomb> {
    static constexpr auto value = libint2::Operator::coulomb;
};

template<>
struct LibintOp<simde::type::el_kinetic> {
    static constexpr auto value = libint2::Operator::kinetic;
};

template<>
struct LibintOp<simde::type::el_identity> {
    static constexpr auto value = libint2::Operator::overlap;
};

template<>
struct LibintOp<simde::type::el_nuc_coulomb> {
    static constexpr auto value = libint2::Operator::nuclear;
};

template<>
struct LibintOp<simde::type::el_el_stg> {
    static constexpr auto value = libint2::Operator::stg;
};

template<>
struct LibintOp<simde::type::el_el_f12_commutator> {
    static constexpr auto value = libint2::Operator::stg;
};

template<>
struct LibintOp<simde::type::el_el_yukawa> {
    static constexpr auto value = libint2::Operator::yukawa;
};

template<>
struct LibintOp<simde::type::el_dipole> {
    static constexpr auto value = libint2::Operator::emultipole1;
};

template<>
struct LibintOp<simde::type::el_quadrupole> {
    static constexpr auto value = libint2::Operator::emultipole2;
};

template<>
struct LibintOp<simde::type::el_octupole> {
    static constexpr auto value = libint2::Operator::emultipole3;
};

} // namespace detail_

template<typename T>
static constexpr auto op_v = detail_::LibintOp<T>::value;

} // namespace integrals::libint
