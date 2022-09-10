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
#include "type_traits.hpp"
#include <libint2.hpp>
#include <simde/types.hpp>

namespace integrals::detail_ {

/** @brief Constructs a Libint engine.
 *
 *  @tparam OpType The type of the integral operator
 *  @param[in] bases A vector of Libint basis sets
 *  @param[in] op A NWX operator for the integral
 *  @param[in] thresh The precision threshold of the integrals
 *  @returns An engine to compute the values of the associated integrals
 */
template<typename OpType>
auto make_engine(const std::vector<libint2::BasisSet>& bases, const OpType& op,
                 double thresh) {
    /// Variables for engine construction
    constexpr auto libint_op = integrals::op_v<OpType>;
    auto max_nprims          = libint2::max_nprim(bases[0]);
    auto max_l               = libint2::max_l(bases[0]);
    std::size_t deriv        = 0;

    /// Find max_nprims and max_l in bases
    for(auto set : bases) {
        max_nprims = std::max(max_nprims, libint2::max_nprim(set));
        max_l      = std::max(max_l, libint2::max_l(set));
    }

    /// Construct engine and handl specialized settings
    if(!libint2::initialized()) libint2::initialize();
    libint2::Engine engine(libint_op, max_nprims, max_l, deriv, thresh);
    /// Libint not acknowleding max_nprim in ctor?
    engine.set_max_nprim(max_nprims);

    if(libint2::rank(libint_op) == 2) {
        if(bases.size() == 2) engine.set(libint2::BraKet::xs_xs);
        if(bases.size() == 3) engine.set(libint2::BraKet::xs_xx);
    }

    if constexpr(std::is_same_v<OpType, simde::type::el_nuc_coulomb>) {
        const auto& nuclei = op.template at<1>();
        std::vector<std::pair<double, std::array<double, 3>>> qs;
        for(const auto& ai : nuclei)
            qs.emplace_back(static_cast<const double&>(ai.Z()), ai.coords());
        engine.set_params(qs);
    } else if constexpr(std::is_same_v<OpType, simde::type::el_el_stg> ||
                        std::is_same_v<OpType, simde::type::el_el_yukawa>) {
        const auto& stg = op.template at<0>();
        engine.set_params(stg.exponent);
    } else if constexpr(std::is_same_v<OpType,
                                       simde::type::el_el_f12_commutator>) {
        const auto& stg = op.template at<0>();
        engine.set_params(2.0 * stg.exponent);
    } else if constexpr(std::is_same_v<OpType, simde::type::el_dipole> ||
                        std::is_same_v<OpType, simde::type::el_quadrupole> ||
                        std::is_same_v<OpType, simde::type::el_octupole>) {
        std::array<double, 3> origin{0, 0, 0};
        engine.set_params(origin);
    }

    return engine;
}

/** @brief Wrap the call of LibInt2 engine so it can take a variable number
 * of shell inputs.
 *
 * @tparam Is A variadic parameter pack of integers from [0,NBases) to
 * expand.
 * @param engine The LibInt2 engine that computes integrals
 * @param bases The bases sets that hold the shells
 * @param shells The index of the requested shell block
 */
template<std::size_t... Is>
void run_engine_(libint2::Engine& engine, std::vector<libint2::BasisSet>& bases,
                 std::vector<std::size_t>& shells, std::index_sequence<Is...>) {
    engine.compute(bases[Is][shells[Is]]...);
}

} // namespace integrals::detail_