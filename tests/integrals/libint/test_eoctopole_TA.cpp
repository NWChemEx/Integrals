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

#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <simde/tensor_representation/tensor_representation.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

using namespace integrals;
using namespace mokup;

TEST_CASE("Octupole") {
    using i_op   = simde::type::el_identity;
    using d_op   = simde::type::el_dipole;
    using q_op   = simde::type::el_quadrupole;
    using o_op   = simde::type::el_octupole;
    using s_type = simde::AOTensorRepresentation<2, i_op>;
    using d_type = simde::AOTensorRepresentation<2, d_op>;
    using q_type = simde::AOTensorRepresentation<2, q_op>;
    using o_type = simde::AOTensorRepresentation<2, o_op>;

    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    const auto name = molecule::h2o;
    const auto bs   = basis_set::sto3g;
    auto aos        = get_bases(name, bs);
    std::vector bases{bs, bs};
    d_op r;
    q_op r2;
    o_op r3;

    /// TODO: this needs to actually test something.
    // SECTION("overlap matrix") {
    //     mm.change_input("EDipole", "Origin", origin);
    //     auto [S] = mm.at("EDipole").run_as<s_type>(aos, aos);
    //     REQUIRE(tensorwrapper::ta_helpers::allclose(S, X));
    // }

    SECTION("dipole matrix") {
        auto [D]  = mm.at("EOctupole").run_as<d_type>(aos, r, aos);
        auto corr = get_ao_data(name, bases, property::dipole);
        REQUIRE(tensorwrapper::tensor::allclose(D, corr));
    }

    SECTION("Quadrupole") {
        auto [Q]  = mm.at("EOctupole").run_as<q_type>(aos, r2, aos);
        auto corr = get_ao_data(name, bases, property::quadrupole);
        REQUIRE(tensorwrapper::tensor::allclose(Q, corr));
    }

    SECTION("Octupole") {
        auto [O]  = mm.at("EOctupole").run_as<o_type>(aos, r3, aos);
        auto corr = get_ao_data(name, bases, property::octopole);
        REQUIRE(tensorwrapper::tensor::allclose(O, corr));
    }
}
