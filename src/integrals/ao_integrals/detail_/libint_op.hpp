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
#include <libint2.hpp>
#include <simde/types.hpp>

namespace integrals::ao_integrals::detail_ {

template<typename T>
struct LibintOp;

template<>
struct LibintOp<simde::type::v_ee_type> {
    static constexpr auto value = libint2::Operator::coulomb;
};

template<>
struct LibintOp<simde::type::t_e_type> {
    static constexpr auto value = libint2::Operator::kinetic;
};

template<>
struct LibintOp<simde::type::v_en_type> {
    static constexpr auto value = libint2::Operator::nuclear;
};

template<>
struct LibintOp<simde::type::s_e_type> {
    static constexpr auto value = libint2::Operator::overlap;
};

template<typename T>
static constexpr auto op_v = detail_::LibintOp<T>::value;

} // namespace integrals::ao_integrals::detail_
