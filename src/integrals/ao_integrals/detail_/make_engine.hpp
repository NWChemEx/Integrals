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
#include "libint_op.hpp"
#include <libint2.hpp>
#include <simde/types.hpp>

namespace integrals::ao_integrals::detail_ {

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
                 double thresh, std::size_t deriv = 0) {
    /// Variables for engine construction
    constexpr auto libint_op = integrals::ao_integrals::detail_::op_v<OpType>;
    auto max_nprims          = libint2::max_nprim(bases[0]);
    auto max_l               = libint2::max_l(bases[0]);

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

    if constexpr(std::is_same_v<OpType, simde::type::v_en_type>) {
        const auto& nuclei = op.rhs_particle();
        std::vector<std::pair<double, std::array<double, 3>>> qs;
        for(const auto& ai : nuclei) {
            std::array<double, 3> coords{ai.x(), ai.y(), ai.z()};
            qs.emplace_back(static_cast<const double&>(ai.Z()), coords);
        }
        engine.set_params(qs);
    }

    return engine;
}

} // namespace integrals::ao_integrals::detail_
